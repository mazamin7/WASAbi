#include "dct_volume.h"
#include <cufft.h>
#include <cuda_runtime.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// -----------------------------------------------------------------------
// DCT-II via 2N even extension + R2C FFT, matching FFTW REDFT10
// -----------------------------------------------------------------------
// Even extension kernel: mirrors original [0..N-1] to [0..2N-1]
// by reflecting around the N-0.5 boundary (half-sample symmetry)
__global__ void EvenExtensionKernel(
    const double* __restrict__ d_in,
    double* __restrict__ d_extended_out,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    int w2 = 2 * w;
    int h2 = 2 * h;
    int d2 = 2 * d;

    if (x < w2 && y < h2 && z < d2) {
        int orig_x = (x < w) ? x : (w2 - 1 - x);
        int orig_y = (y < h) ? y : (h2 - 1 - y);
        int orig_z = (z < d) ? z : (d2 - 1 - z);
        double val = d_in[orig_z * h * w + orig_y * w + orig_x];
        d_extended_out[z * h2 * w2 + y * w2 + x] = val;
    }
}

// Post-process DCT-II: extract DCT modes from R2C FFT of even-extended signal.
// The R2C transform of the 2N even extended signal gives:
//   Y[k] = FFT[ext][k]
// The DCT-II modes are:
//   X_dct[k] = Re( Y[k] * exp(-j*pi*k/(2N)) ) / (2*N) * (unnormalized fftw convention)
// 
// Working backwards from FFTW REDFT10 which scales by 2N and we then divide by
// 2*sqrt(2*N) for orthonormal: the full internal scale is 2*sqrt(2*N).
// So here: mode[k] = Re(Y[k] * twiddle[k]) / (2*sqrt(2*N_x*N_y*N_z))
__global__ void PostProcessDctIIKernel(
    const cufftDoubleComplex* __restrict__ d_complex_in,
    double* __restrict__ d_modes_out,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int w2 = 2 * w;
        int h2 = 2 * h;
        int d2 = 2 * d;
        // R2C output: layout is (d2, h2, w2/2+1) = (d2, h2, w+1)
        int Nhalf_x = w + 1;
        int idx_complex = z * h2 * Nhalf_x + y * Nhalf_x + x;

        cufftDoubleComplex G = d_complex_in[idx_complex];

        // Twiddle per dimension: exp(-j*pi*k/(2*N))
        // 3D separable: twiddle = twiddle_x * twiddle_y * twiddle_z
        // = exp(-j * pi/2 * (x/w + y/h + z/d))
        double theta_x = -M_PI * x / (2.0 * w);
        double theta_y = -M_PI * y / (2.0 * h);
        double theta_z = -M_PI * z / (2.0 * d);
        double theta = theta_x + theta_y + theta_z;

        double con = cos(theta);
        double son = sin(theta);

        // Re( (Gx + j*Gy) * (cos + j*sin) ) = Gx*cos - Gy*sin
        double dct_coeff = G.x * con - G.y * son;

        // Match WASAbi DctVolume.cpp Normalization:
        // CPU code does: modes[idx] = fftw_output[idx] / (2.0 * sqrt(2.0 * N_total))
        // Our G is the result of R2C FFT of 2N even extension.
        // For REDFT10: X_dct[k] = Re( FFT_2N[k] * exp(-j*pi*k/(2N)) )
        // BUT the CPU uses REDFT10 which scales by 2N *internally* if not using FFTW_MEASURE? 
        // No, standard FFTW r2r REDFT10 is unnormalized.
        // Let's use the exact scaling: coeff = G * twiddle, then divide by 2*sqrt(2*N)
        // Wait, the R2C on 2N signal also has a 2N scale factor? No, cuFFT R2C is unnormalized too.
        // Standard D2Z of length M signal has a sum of real values.
        // Total normalization should bring us to the same magnitude as CPU.
        double scale = 2.0 * sqrt(2.0 * (double)(w * h * d));
        d_modes_out[z * h * w + y * w + x] = dct_coeff / scale;
    }
}

// Pre-process for IDCT (DCT-III): multiply modes by conjugate twiddle and insert into
// the Hermitian-symmetric complex array for C2R.
// We fill the first octant (x<w, y<h, z<d) with modes * conj_twiddle.
// The C2R transform expects Hermitian symmetry which we enforce by only filling
// the half-spectrum (x <= w) and relying on cuFFT's C2R convention.
__global__ void PreProcessIdctKernel(
    const double* __restrict__ d_modes_in,
    cufftDoubleComplex* __restrict__ d_complex_out,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    int w2 = 2 * w;
    int h2 = 2 * h;
    int d2 = 2 * d;
    int Nhalf_x = w + 1;

    if (x < Nhalf_x && y < h2 && z < d2) {
        cufftDoubleComplex G;
        G.x = 0.0;
        G.y = 0.0;

        int folded_y = (y <= h) ? y : (h2 - y);
        int folded_z = (z <= d) ? z : (d2 - z);
        
        if (x < w && folded_y < h && folded_z < d) {
            double mode_val = d_modes_in[folded_z * h * w + folded_y * w + x];
            
            double sign_y = (y < h) ? 1.0 : -1.0;
            double sign_z = (z < d) ? 1.0 : -1.0;
            
            double theta_x = +M_PI * x / (2.0 * w); // Always positive (x < w)
            double theta_y = sign_y * M_PI * folded_y / (2.0 * h);
            double theta_z = sign_z * M_PI * folded_z / (2.0 * d);
            
            double theta = theta_x + theta_y + theta_z;
            
            G.x = mode_val * cos(theta);
            G.y = mode_val * sin(theta);
        }

        d_complex_out[z * h2 * Nhalf_x + y * Nhalf_x + x] = G;
    }
}

// Post-process after C2R IDCT: extract first N values and divide by 2N
__global__ void PostProcessIdctKernel(
    const double* __restrict__ d_c2r_out,
    double* __restrict__ d_values_out,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int w2 = 2 * w;
        int h2 = 2 * h;
        // The C2R DOES NOT divide by anything in cuFFT.
        double scale = 2.0 * sqrt(2.0 * (double)(w * h * d));
        double val = d_c2r_out[z * h2 * w2 + y * w2 + x];
        d_values_out[z * h * w + y * w + x] = val / scale; 
    }
}

// -----------------------------------------------------------------------
// DctVolume implementation
// -----------------------------------------------------------------------
DctVolume::DctVolume(int w, int h, int d)
    : width_(w), height_(h), depth_(d)
    , d_values_(nullptr), d_modes_(nullptr), d_complex_modes_(nullptr)
{
    size_t size_orig = (size_t)width_ * height_ * depth_ * sizeof(double);
    cudaMalloc((void**)&d_values_, size_orig);
    cudaMalloc((void**)&d_modes_,  size_orig);
    cudaMemset(d_values_, 0, size_orig);
    cudaMemset(d_modes_,  0, size_orig);

    int w2 = 2 * width_;
    int h2 = 2 * height_;
    int d2 = 2 * depth_;

    // R2C output: (d2, h2, w2/2+1)
    size_t size_complex = (size_t)d2 * h2 * (width_ + 1) * sizeof(cufftDoubleComplex);
    cudaMalloc((void**)&d_complex_modes_, size_complex);

    // Extended real buffer for DCT
    size_t size_ext = (size_t)w2 * h2 * d2 * sizeof(double);
    cudaMalloc((void**)&d_extended_, size_ext);

    // D2Z = real-to-complex, Z2D = complex-to-real
    cufftPlan3d(&r2c_plan_, d2, h2, w2, CUFFT_D2Z);
    cufftPlan3d(&c2r_plan_, d2, h2, w2, CUFFT_Z2D);
}

DctVolume::~DctVolume()
{
    cufftDestroy(r2c_plan_);
    cufftDestroy(c2r_plan_);
    cudaFree(d_values_);
    cudaFree(d_modes_);
    cudaFree(d_complex_modes_);
    cudaFree(d_extended_);
}

void DctVolume::ExecuteDct()
{
    int w2 = 2 * width_;
    int h2 = 2 * height_;
    int d2 = 2 * depth_;

    dim3 blockSize(8, 8, 8);
    dim3 extGridSize((w2 + blockSize.x - 1) / blockSize.x,
                     (h2 + blockSize.y - 1) / blockSize.y,
                     (d2 + blockSize.z - 1) / blockSize.z);

    // 1. Even extension
    EvenExtensionKernel<<<extGridSize, blockSize>>>(d_values_, d_extended_, width_, height_, depth_);

    // 2. R2C FFT on the 2N extended signal
    cufftExecD2Z(r2c_plan_, d_extended_, d_complex_modes_);

    // 3. Post-process to extract DCT-II modes
    dim3 origGridSize((width_  + blockSize.x - 1) / blockSize.x,
                      (height_ + blockSize.y - 1) / blockSize.y,
                      (depth_  + blockSize.z - 1) / blockSize.z);
    PostProcessDctIIKernel<<<origGridSize, blockSize>>>(d_complex_modes_, d_modes_, width_, height_, depth_);
}

void DctVolume::ExecuteIdct()
{
    int w2 = 2 * width_;
    int h2 = 2 * height_;
    int d2 = 2 * depth_;

    // 1. Pre-process: fill Hermitian complex input from modes
    size_t size_complex = (size_t)d2 * h2 * (width_ + 1) * sizeof(cufftDoubleComplex);
    cudaMemset(d_complex_modes_, 0, size_complex);

    dim3 blockSize(8, 8, 8);
    int Nhalf_x = width_ + 1;
    dim3 preGridSize((Nhalf_x   + blockSize.x - 1) / blockSize.x,
                     (h2        + blockSize.y - 1) / blockSize.y,
                     (d2        + blockSize.z - 1) / blockSize.z);
    PreProcessIdctKernel<<<preGridSize, blockSize>>>(d_modes_, d_complex_modes_, width_, height_, depth_);

    // 2. C2R FFT
    cufftExecZ2D(c2r_plan_, d_complex_modes_, d_extended_);

    // 3. Extract and normalize
    dim3 origGridSize((width_  + blockSize.x - 1) / blockSize.x,
                      (height_ + blockSize.y - 1) / blockSize.y,
                      (depth_  + blockSize.z - 1) / blockSize.z);
    PostProcessIdctKernel<<<origGridSize, blockSize>>>(d_extended_, d_values_, width_, height_, depth_);
}

void DctVolume::reset()
{
    size_t size = (size_t)width_ * height_ * depth_ * sizeof(double);
    cudaMemset((void*)d_values_, 0, size);
    cudaMemset((void*)d_modes_,  0, size);
}
