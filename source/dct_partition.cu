#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// The massively parallel O(N^3) physics update kernel - now optimized with pre-calculated coefficients
__global__ void DctUpdateKernel(
    double* __restrict__ d_pressure_modes,
    double* __restrict__ d_velocity_modes,
    const double* __restrict__ d_force_modes,
    const double* __restrict__ d_A,
    const double* __restrict__ d_B,
    const double* __restrict__ d_C,
    const double* __restrict__ d_D,
    const double* __restrict__ d_E,
    const double* __restrict__ d_F,
    int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int idx = z * h * w + y * w + x;
        
        double pres = d_pressure_modes[idx];
        double vel = d_velocity_modes[idx];
        double force = d_force_modes[idx];
        
        // Optimized physics update using pre-calculated coefficients
        // next_vel = A*vel + B*pres + C*force
        // next_pres = D*vel + E*pres + F*force
        double next_vel  = d_A[idx] * vel + d_B[idx] * pres + d_C[idx] * force;
        double next_pres = d_D[idx] * vel + d_E[idx] * pres + d_F[idx] * force;

        d_velocity_modes[idx] = next_vel;
        d_pressure_modes[idx] = next_pres;
    }
}

// Kernel to pre-calculate w0, alpha, AND physics coefficients natively on the GPU during initialization
__global__ void InitConstantsKernel(
    double* __restrict__ d_w0, double* __restrict__ d_alpha,
    double* __restrict__ d_A, double* __restrict__ d_B, double* __restrict__ d_C,
    double* __restrict__ d_D, double* __restrict__ d_E, double* __restrict__ d_F,
    int w, int h, int d,
    double lx2, double ly2, double lz2,
    double c0, double a1, double a2, double dt)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x; // width (x)
    int j = blockIdx.y * blockDim.y + threadIdx.y; // height (y)
    int i = blockIdx.z * blockDim.z + threadIdx.z; // depth (z)

    if (k < w && j < h && i < d) {
        int idx = i * h * w + j * w + k;
        double w0 = c0 * M_PI * sqrt(i * i / lz2 + j * j / ly2 + k * k / lx2);
        double alpha = a1 + a2 * w0 * w0;
        
        d_w0[idx] = w0;
        d_alpha[idx] = alpha;

        double A = 1.0, B = 0.0, C = 0.0;
        double D = 0.0, E = 1.0, F = 0.0;

        if ((idx == 0) && (alpha == 0.0)) {
            A = 1.0; B = 0.0; C = dt;
            D = dt;  E = 1.0; F = (dt * dt / 2.0);
        }
        else if ((idx == 0) && (alpha > 0.0)) {
            double e2at = exp(-2.0 * alpha * dt);
            A = e2at; B = 0.0; C = (1.0 - e2at) / (2.0 * alpha);
            D = (1.0 - e2at) / (2.0 * alpha); E = 1.0; F = ((e2at - 1.0) / (4.0 * alpha * alpha) + 1.0 / (2.0 * alpha) * dt);
        }
        else if ((idx > 0) && (alpha < w0)) {
            double inv_w02 = 1.0 / (w0 * w0);
            double alpha_sqr = alpha * alpha;
            double omega = sqrt(w0 * w0 - alpha_sqr);
            double cwt = cos(omega * dt);
            double swt = sin(omega * dt);
            double eatm = exp(-alpha * dt);
            double inv_w = 1.0 / omega;

            A = eatm * (cwt - alpha * inv_w * swt);
            B = -eatm * (omega + alpha_sqr * inv_w) * swt;
            C = -B * inv_w02;

            D = eatm * swt * inv_w;
            E = eatm * (cwt + alpha * inv_w * swt);
            F = (1.0 - E) * inv_w02;
        }
        else if ((idx > 0) && (alpha > w0)) {
            double inv_w02 = 1.0 / (w0 * w0);
            double alpha_sqr = alpha * alpha;
            double alphad = sqrt(alpha_sqr - w0 * w0);
            double alpha1 = alpha + alphad;
            double alpha2 = alpha - alphad;
            double eat1 = exp(-alpha1 * dt);
            double eat2 = exp(-alpha2 * dt);

            A = (0.5 * (eat1 + eat2) + 0.5 / alphad * alpha * (eat1 - eat2));
            B = (-0.5 * (alpha1 * eat1 + alpha2 * eat2) - 0.5 / alphad * alpha * (alpha2 * eat2 - alpha1 * eat1));
            C = inv_w02 * 0.5 * (alpha1 * eat1 + alpha2 * eat2 + alpha / alphad * (alpha2 * eat2 - alpha1 * eat1));

            D = 0.5 / alphad * (eat2 - eat1);
            E = (eat1 + eat2 + 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1));
            F = inv_w02 * (1.0 - (eat1 + eat2 + 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1)));
        }

        d_A[idx] = A; d_B[idx] = B; d_C[idx] = C;
        d_D[idx] = D; d_E[idx] = E; d_F[idx] = F;
    }
}

DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d)
    : Partition(xs, ys, zs, w, h, d)
{
    // The sizes
    size_t vol_size = w * h * d * sizeof(double);
    cudaMalloc((void**)&d_w0_, vol_size);
    cudaMalloc((void**)&d_alpha_, vol_size);
    
    cudaMalloc((void**)&d_coef_A_, vol_size);
    cudaMalloc((void**)&d_coef_B_, vol_size);
    cudaMalloc((void**)&d_coef_C_, vol_size);
    cudaMalloc((void**)&d_coef_D_, vol_size);
    cudaMalloc((void**)&d_coef_E_, vol_size);
    cudaMalloc((void**)&d_coef_F_, vol_size);

    // Swap the class members to use the DctVolume
    pressure_vol_ = new DctVolume(w, h, d);
    velocity_vol_ = new DctVolume(w, h, d);
    force_vol_ = new DctVolume(w, h, d);
    
    // Wire the base Partition pointers to point to the Volume device arrays
    d_pressure_ = pressure_vol_->d_values_;
    d_velocity_ = velocity_vol_->d_values_;
    d_force_ = force_vol_->d_values_;

    InitializeConstants();
}

DctPartition::~DctPartition()
{
    cudaFree(d_w0_);
    cudaFree(d_alpha_);
    cudaFree(d_coef_A_); cudaFree(d_coef_B_); cudaFree(d_coef_C_);
    cudaFree(d_coef_D_); cudaFree(d_coef_E_); cudaFree(d_coef_F_);
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
        d_coef_A_, d_coef_B_, d_coef_C_, 
        d_coef_D_, d_coef_E_, d_coef_F_, 
        width_, height_, depth_, lx2, ly2, lz2, c0_, air_absorption_alpha1_, air_absorption_alpha2_, dt_);
}

void DctPartition::Update()
{
    // 1. Transform space to frequency (DCT)
    pressure_vol_->ExecuteDct(stream_);
    velocity_vol_->ExecuteDct(stream_);
    force_vol_->ExecuteDct(stream_);

    // 2. Optimized Physics Update using pre-calculated coefficients
    dim3 blockSize(8, 8, 8);
    dim3 gridSize((width_ + blockSize.x - 1) / blockSize.x,
                  (height_ + blockSize.y - 1) / blockSize.y,
                  (depth_ + blockSize.z - 1) / blockSize.z);

    DctUpdateKernel<<<gridSize, blockSize, 0, stream_>>>(
        pressure_vol_->d_modes_,
        velocity_vol_->d_modes_,
        force_vol_->d_modes_,
        d_coef_A_, d_coef_B_, d_coef_C_,
        d_coef_D_, d_coef_E_, d_coef_F_,
        width_, height_, depth_
    );

    // 3. Transform frequency back to space (IDCT)
    velocity_vol_->ExecuteIdct(stream_);
    pressure_vol_->ExecuteIdct(stream_);
}

