#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// The massively parallel O(N^3) physics update kernel
__global__ void DctUpdateKernel(
    double* __restrict__ d_pressure_modes,
    double* __restrict__ d_velocity_modes,
    const double* __restrict__ d_force_modes,
    const double* __restrict__ d_w0,
    const double* __restrict__ d_alpha,
    double dt, int w, int h, int d)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < w && y < h && z < d) {
        int idx = z * h * w + y * w + x;
        double w0 = d_w0[idx];
        double alpha = d_alpha[idx];
        
        double pres = d_pressure_modes[idx];
        double vel = d_velocity_modes[idx];
        double force = d_force_modes[idx];
        
        double next_vel, next_pres;

        if ((idx == 0) && (alpha == 0.0)) {
            next_vel = vel + dt * force;
            next_pres = dt * vel + pres + (dt * dt / 2.0) * force;
        }
        else if ((idx == 0) && (alpha > 0.0)) {
            double e2at = exp(-2.0 * alpha * dt);
            next_vel = e2at * vel + (1.0 - e2at) / (2.0 * alpha) * force;
            next_pres = (1.0 - e2at) / (2.0 * alpha) * vel + pres + ((e2at - 1.0) / (4.0 * alpha * alpha) + 1.0 / (2.0 * alpha) * dt) * force;
        }
        else if ((idx > 0) && (alpha < w0)) {
            double inv_w02 = 1.0 / (w0 * w0);
            double alpha_sqr = alpha * alpha;
            double w = sqrt(w0 * w0 - alpha_sqr);
            double cwt = cos(w * dt);
            double swt = sin(w * dt);
            double eatm = exp(-alpha * dt);
            double inv_w = 1.0 / w;
            double xe = force * inv_w02;

            next_vel = eatm * (vel * (cwt - alpha * inv_w * swt) - (w + alpha_sqr * inv_w) * (pres - xe) * swt);
            next_pres = xe + eatm * ((pres - xe) * (cwt + alpha * inv_w * swt) + swt * inv_w * vel);
        }
        else if ((idx > 0) && (alpha > w0)) {
            double inv_w02 = 1.0 / (w0 * w0);
            double alpha_sqr = alpha * alpha;
            double alphad = sqrt(alpha_sqr - w0 * w0);
            double alpha1 = alpha + alphad;
            double alpha2 = alpha - alphad;
            double eat1 = exp(-alpha1 * dt);
            double eat2 = exp(-alpha2 * dt);

            next_vel = (0.5 * (eat1 + eat2) + 0.5 / alphad * alpha * (eat1 - eat2)) * vel 
                         + (-0.5 * (alpha1 * eat1 + alpha2 * eat2) - 0.5 / alphad * alpha * (alpha2 * eat2 - alpha1 * eat1)) * pres 
                         + inv_w02 * 0.5 * (alpha1 * eat1 + alpha2 * eat2 + alpha / alphad * (alpha2 * eat2 - alpha1 * eat1)) * force;
            
            next_pres = 0.5 / alphad * (eat2 - eat1) * vel 
                          + (eat1 + eat2 + 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1)) * pres 
                          + inv_w02 * (1.0 - alpha1 - alpha2 - 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1)) * force;
        }
        else {
            next_vel = vel; 
            next_pres = pres;
        }

        d_velocity_modes[idx] = next_vel;
        d_pressure_modes[idx] = next_pres;
    }
}

// Kernel to pre-calculate w0 and alpha natively on the GPU during initialization
__global__ void InitConstantsKernel(
    double* __restrict__ d_w0,
    double* __restrict__ d_alpha,
    int w, int h, int d,
    double lx2, double ly2, double lz2,
    double c0, double a1, double a2)
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
    }
}

DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d)
    : Partition(xs, ys, zs, w, h, d)
{
    // The sizes
    size_t vol_size = w * h * d * sizeof(double);
    cudaMalloc((void**)&d_w0_, vol_size);
    cudaMalloc((void**)&d_alpha_, vol_size);

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

    InitConstantsKernel<<<gridSize, blockSize>>>(d_w0_, d_alpha_, width_, height_, depth_, lx2, ly2, lz2, c0_, air_absorption_alpha1_, air_absorption_alpha2_);
}

void DctPartition::Update()
{
    // 1. Transform space to frequency (DCT)
    pressure_vol_->ExecuteDct();
    velocity_vol_->ExecuteDct();
    force_vol_->ExecuteDct();

    // 2. O(N^3) Physics Update mapped to millions of CUDA threads
    dim3 blockSize(8, 8, 8);
    dim3 gridSize((width_ + blockSize.x - 1) / blockSize.x,
                  (height_ + blockSize.y - 1) / blockSize.y,
                  (depth_ + blockSize.z - 1) / blockSize.z);

    DctUpdateKernel<<<gridSize, blockSize>>>(
        pressure_vol_->d_modes_,
        velocity_vol_->d_modes_,
        force_vol_->d_modes_,
        d_w0_,
        d_alpha_,
        dt_, width_, height_, depth_
    );

    // 3. Transform frequency back to space (IDCT)
    velocity_vol_->ExecuteIdct();
    pressure_vol_->ExecuteIdct();
}

