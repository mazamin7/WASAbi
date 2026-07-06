#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Kernel to pre-calculate w0, alpha, AND physics coefficients natively on the GPU during initialization
__global__ void InitConstantsKernel(
    double* __restrict__ d_w0, double* __restrict__ d_alpha,
    double* __restrict__ d_S11, double* __restrict__ d_S12, double* __restrict__ d_S21, double* __restrict__ d_S22,
    int w, int h, int d,
    double lx2, double ly2, double lz2,
    double c0, double a1, double a2, double dt)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x; // width (x)
    int j = blockIdx.y * blockDim.y + threadIdx.y; // height (y)
    int i = blockIdx.z * blockDim.z + threadIdx.z; // depth (z)

    if (k < w && j < h && i < d) {
        int idx = i * h * w + j * w + k;
        double w0 = c0 * M_PI * sqrt((double)(i * i) / lz2 + (double)(j * j) / ly2 + (double)(k * k) / lx2);
        double alpha = a1 + a2 * w0 * w0;
        
        d_w0[idx] = w0;
        d_alpha[idx] = alpha;

        double w02 = w0 * w0;
        double S11 = 1.0, S12 = dt, S21 = 0.0, S22 = 1.0;

        if (idx == 0) {
            // Constant mode (\lambda_m = 0 -> w0 = 0)
            if (alpha == 0.0) {
                S11 = 1.0; S12 = dt;
                S21 = 0.0; S22 = 1.0;
            } else {
                double e_alpha_dt = exp(-alpha * dt);
                S11 = 1.0; S12 = (1.0 - e_alpha_dt) / alpha;
                S21 = 0.0; S22 = e_alpha_dt;
            }
        } else {
            // Non-constant modes (\lambda_m < 0 -> w0 > 0)
            double delta_sq = alpha * alpha - 4.0 * w02;
            double E = exp(-alpha * dt / 2.0);
            double Ch = 0.0, Sh = 0.0;

            if (delta_sq < 0.0) {
                // Underdamped
                double omega = sqrt(-delta_sq) / 2.0;
                Ch = cos(omega * dt);
                Sh = sin(omega * dt) / omega;
            } else if (delta_sq == 0.0) {
                // Critically damped
                Ch = 1.0;
                Sh = dt;
            } else {
                // Overdamped
                double delta = sqrt(delta_sq);
                Ch = cosh(delta * dt / 2.0);
                Sh = sinh(delta * dt / 2.0) / (delta / 2.0);
            }

            S11 = E * (Ch + Sh * alpha / 2.0);
            S12 = E * Sh;
            S21 = E * Sh * (-w02);
            S22 = E * (Ch - Sh * alpha / 2.0);
        }

        d_S11[idx] = S11; d_S12[idx] = S12; 
        d_S21[idx] = S21; d_S22[idx] = S22;
    }
}

__global__ void DctDriftKernel(
    double* __restrict__ d_pressure_modes,
    double* __restrict__ d_velocity_modes,
    const double* __restrict__ d_S11,
    const double* __restrict__ d_S12,
    const double* __restrict__ d_S21,
    const double* __restrict__ d_S22,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int idx = z * h * w + y * w + x;
        
        double pres = d_pressure_modes[idx];
        double vel = d_velocity_modes[idx];
        
        // Exact drift
        double next_pres = d_S11[idx] * pres + d_S12[idx] * vel;
        double next_vel  = d_S21[idx] * pres + d_S22[idx] * vel;

        d_pressure_modes[idx] = next_pres;
        d_velocity_modes[idx] = next_vel;
    }
}

__global__ void DctKickKernel(
    double* __restrict__ d_velocity_modes,
    const double* __restrict__ d_force_modes,
    double half_dt,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int idx = z * h * w + y * w + x;
        d_velocity_modes[idx] += half_dt * d_force_modes[idx];
    }
}

__global__ void MergeResiduesIntoForceKernel(double* d_force, const double* d_residue, int size) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) d_force[i] += d_residue[i];
}


DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d)
    : Partition(xs, ys, zs, w, h, d, false) // Don't allocate default Partition buffers (prevents leak)
{
    // The sizes
    size_t vol_size = (size_t)w * h * d * sizeof(double);
    cudaMalloc((void**)&d_w0_, vol_size);
    cudaMalloc((void**)&d_alpha_, vol_size);
    
    cudaMalloc((void**)&d_S11_, vol_size);
    cudaMalloc((void**)&d_S12_, vol_size);
    cudaMalloc((void**)&d_S21_, vol_size);
    cudaMalloc((void**)&d_S22_, vol_size);

    // Allocate shared buffers for DCT operations to save VRAM
    int w2 = 2 * w, h2 = 2 * h, d2 = 2 * d;
    size_t size_ext = (size_t)w2 * h2 * d2 * sizeof(double);
    size_t size_complex = (size_t)d2 * h2 * (w + 1) * sizeof(cufftDoubleComplex);
    
    cudaMalloc((void**)&d_shared_ext_, size_ext);
    cudaMalloc((void**)&d_shared_complex_, size_complex);

    // Create DctVolumes using the shared buffers
    pressure_vol_ = new DctVolume(w, h, d, d_shared_ext_, d_shared_complex_);
    velocity_vol_ = new DctVolume(w, h, d, d_shared_ext_, d_shared_complex_);
    force_vol_    = new DctVolume(w, h, d, d_shared_ext_, d_shared_complex_);
    
    // Wire the base Partition pointers to point to the Volume device arrays
    d_pressure_ = pressure_vol_->d_values_;
    d_velocity_ = velocity_vol_->d_values_;
    d_force_    = force_vol_->d_values_;

    InitializeConstants();
}

DctPartition::~DctPartition()
{
    cudaFree(d_shared_ext_);
    cudaFree(d_shared_complex_);
    cudaFree(d_w0_);
    cudaFree(d_alpha_);
    cudaFree(d_S11_); cudaFree(d_S12_); cudaFree(d_S21_); cudaFree(d_S22_);
    delete pressure_vol_;
    delete velocity_vol_;
    delete force_vol_;
}

void DctPartition::InitializeConstants()
{
    double lx2 = width_ * width_ * dh_ * dh_;
    double ly2 = height_ * height_ * dh_ * dh_;
    double lz2 = depth_ * depth_ * dh_ * dh_;
    dim3 blockSize(8, 8, 8);
    dim3 gridSize((width_ + blockSize.x - 1) / blockSize.x,
                  (height_ + blockSize.y - 1) / blockSize.y,
                  (depth_ + blockSize.z - 1) / blockSize.z);

    InitConstantsKernel<<<gridSize, blockSize>>>(
        d_w0_, d_alpha_, 
        d_S11_, d_S12_, d_S21_, d_S22_, 
        width_, height_, depth_, lx2, ly2, lz2, c0_, air_absorption_alpha1_, air_absorption_alpha2_, dt_);
}

void DctPartition::HalfKick()
{
    dim3 blockSize(8, 8, 8);
    dim3 gridSize((width_ + blockSize.x - 1) / blockSize.x,
                  (height_ + blockSize.y - 1) / blockSize.y,
                  (depth_ + blockSize.z - 1) / blockSize.z);

    DctKickKernel<<<gridSize, blockSize, 0, stream_>>>(
        velocity_vol_->d_modes_,
        force_vol_->d_modes_,
        dt_ / 2.0,
        width_, height_, depth_
    );
}

void DctPartition::Drift()
{
    dim3 blockSize(8, 8, 8);
    dim3 gridSize((width_ + blockSize.x - 1) / blockSize.x,
                  (height_ + blockSize.y - 1) / blockSize.y,
                  (depth_ + blockSize.z - 1) / blockSize.z);

    DctDriftKernel<<<gridSize, blockSize, 0, stream_>>>(
        pressure_vol_->d_modes_,
        velocity_vol_->d_modes_,
        d_S11_, d_S12_, d_S21_, d_S22_,
        width_, height_, depth_
    );
}

void DctPartition::MergeResiduesIntoForce()
{
    int size = width_ * height_ * depth_;
    int threads = 256;
    int blocks = (size + threads - 1) / threads;
    MergeResiduesIntoForceKernel<<<blocks, threads, 0, stream_>>>(d_force_, d_residue_, size);
}

void DctPartition::Update()
{
    // Deprecated for Strang splitting
}
