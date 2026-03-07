#include "pml_partition.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <cmath>
#include "simulation.h"

// Macro for CUDA error checking
#define CHECK_CUDA(call) \
    { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA Error in %s at line %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    }

__device__ int GetPmlIndex(int x, int y, int z, int width, int height, int depth)
{
    if (x < 0) x = -x - 1;
    else if (x >= width) x = width - 1 - (x - width);
    if (y < 0) y = -y - 1;
    else if (y >= height) y = height - 1 - (y - height);
    if (z < 0) z = -z - 1;
    else if (z >= depth) z = depth - 1 - (z - depth);

    x = (x < 0) ? 0 : ((x >= width) ? width - 1 : x);
    y = (y < 0) ? 0 : ((y >= height) ? height - 1 : y);
    z = (z < 0) ? 0 : ((z >= depth) ? depth - 1 : z);

    return z * height * width + y * width + x;
}

__global__ void PmlUpdateKernel(
    double* p, double* p_new, double* v, double* v_new, 
    double* psi, double* psi_new,
    double* phix, double* phix_new,
    double* phiy, double* phiy_new,
    double* phiz, double* phiz_new,
    const double* zetax, const double* zetay, const double* zetaz,
    const double* force,
    int width, int height, int depth,
    double dh, double dt, double c0)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < width && y < height && z < depth) {
        int idx = GetPmlIndex(x, y, z, width, height, depth);
        
        // 1. Pressure Update
        p_new[idx] = p[idx] + dt * v[idx];
        
        // 2. Auxiliary field psi
        psi_new[idx] = psi[idx] + dt * p_new[idx];

        // 3. Stencil for Delta p
        double coefs[7] = { 2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0 };
        double d2udx2 = 0, d2udy2 = 0, d2udz2 = 0;
        for (int m = 0; m < 7; m++) {
            d2udx2 += coefs[m] * p_new[GetPmlIndex(x + m - 3, y, z, width, height, depth)];
            d2udy2 += coefs[m] * p_new[GetPmlIndex(x, y + m - 3, z, width, height, depth)];
            d2udz2 += coefs[m] * p_new[GetPmlIndex(x, y, z + m - 3, width, height, depth)];
        }
        d2udx2 /= (180.0 * dh * dh);
        d2udy2 /= (180.0 * dh * dh);
        d2udz2 /= (180.0 * dh * dh);

        // 4. Stencil for Nabla phi
        double fourthCoefs[5] = { 1.0, -8.0, 0.0, 8.0, -1.0 };
        double dphidx = 0, dphidy = 0, dphidz = 0;
        for (int m = 0; m < 5; m++) {
            dphidx += fourthCoefs[m] * phix[GetPmlIndex(x + m - 2, y, z, width, height, depth)];
            dphidy += fourthCoefs[m] * phiy[GetPmlIndex(x, y + m - 2, z, width, height, depth)];
            dphidz += fourthCoefs[m] * phiz[GetPmlIndex(x, y, z + m - 2, width, height, depth)];
        }
        dphidx /= (12.0 * dh);
        dphidy /= (12.0 * dh);
        dphidz /= (12.0 * dh);

        // 5. Velocity Update
        // Match PmlPartition.cpp exactly
        double term3 = c0 * c0 * (d2udx2 + d2udy2 + d2udz2);
        double term4 = -(zetax[idx] + zetay[idx] + zetaz[idx]) * v[idx] / 2.0;
        double term5 = -(zetax[idx] * zetay[idx] + zetay[idx] * zetaz[idx] + zetax[idx] * zetaz[idx]) * p_new[idx];
        double term6 = dphidx + dphidy + dphidz;
        // In PmlPartition.cpp line 285: (psi_new + psi)*0.5
        double term7 = -zetax[idx] * zetay[idx] * zetaz[idx] * (psi_new[idx] + psi[idx]) * 0.5;

        v_new[idx] = v[idx] + dt * (term3 + term4 + term5 + term6 + term7 + force[idx]);
        v_new[idx] /= (1.0 + dt / 2.0 * (zetax[idx] + zetay[idx] + zetaz[idx]));

        // 6. Phi Update
        double dudx = 0, dudy = 0, dudz = 0;
        double dpsidx = 0, dpsidy = 0, dpsidz = 0;
        for (int m = 0; m < 5; m++) {
            int ix = GetPmlIndex(x + m - 2, y, z, width, height, depth);
            int iy = GetPmlIndex(x, y + m - 2, z, width, height, depth);
            int iz = GetPmlIndex(x, y, z + m - 2, width, height, depth);
            dudx += fourthCoefs[m] * p_new[ix];
            dudy += fourthCoefs[m] * p_new[iy];
            dudz += fourthCoefs[m] * p_new[iz];
            // Wait, dpsidx in original code was using p_new! 
            // Looking at the code: dpsidx /= (12.0 * dh); using p_new. 
            // In the original PmlPartition.cpp line 316, it uses p_new. 
            // Wait, that looks like a bug in original code (should be psi_new?).
            // Let's stick to original behavior for now.
            dpsidx += fourthCoefs[m] * p_new[ix]; 
            dpsidy += fourthCoefs[m] * p_new[iy];
            dpsidz += fourthCoefs[m] * p_new[iz];
        }
        dudx /= (12.0 * dh); dudy /= (12.0 * dh); dudz /= (12.0 * dh);
        dpsidx /= (12.0 * dh); dpsidy /= (12.0 * dh); dpsidz /= (12.0 * dh);

        phix_new[idx] = phix[idx] + dt * (-zetax[idx] * phix[idx] + c0 * c0 * (zetay[idx] + zetaz[idx] - zetax[idx]) * dudx + c0 * c0 * zetay[idx] * zetaz[idx] * dpsidx);
        phiy_new[idx] = phiy[idx] + dt * (-zetay[idx] * phiy[idx] + c0 * c0 * (zetax[idx] + zetaz[idx] - zetay[idx]) * dudy + c0 * c0 * zetax[idx] * zetaz[idx] * dpsidy);
        phiz_new[idx] = phiz[idx] + dt * (-zetaz[idx] * phiz[idx] + c0 * c0 * (zetax[idx] + zetay[idx] - zetaz[idx]) * dudz + c0 * c0 * zetax[idx] * zetay[idx] * dpsidz);
    }
}

PmlPartition::PmlPartition(std::shared_ptr<Partition> neighbor_part, PmlType type, int xs, int ys, int zs, int w, int h, int d)
    : Partition(xs, ys, zs, w, h, d), type_((PmlPartition::PmlType)type)
{
    info_.type = "PML (CUDA)";
    should_render_ = false;

    if (type_ == P_LEFT || type_ == P_RIGHT) is_x_pml_ = true;
    if (type_ == P_TOP || type_ == P_BOTTOM) is_y_pml_ = true;
    if (type_ == P_FRONT || type_ == P_BACK) is_z_pml_ = true;

    int size = width_ * height_ * depth_;
    size_t bytes = size * sizeof(double);

    CHECK_CUDA(cudaMalloc(&d_psi_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_x_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_y_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_z_, bytes));
    CHECK_CUDA(cudaMalloc(&d_zetax_, bytes));
    CHECK_CUDA(cudaMalloc(&d_zetay_, bytes));
    CHECK_CUDA(cudaMalloc(&d_zetaz_, bytes));
    CHECK_CUDA(cudaMalloc(&d_p_new_, bytes));
    CHECK_CUDA(cudaMalloc(&d_v_new_, bytes));
    CHECK_CUDA(cudaMalloc(&d_psi_new_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_x_new_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_y_new_, bytes));
    CHECK_CUDA(cudaMalloc(&d_phi_z_new_, bytes));

    CHECK_CUDA(cudaMemset(d_psi_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_phi_x_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_phi_y_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_phi_z_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_zetax_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_zetay_, 0, bytes));
    CHECK_CUDA(cudaMemset(d_zetaz_, 0, bytes));

    // Initialize zeta on host and copy to device
    double thickness = Simulation::n_pml_layers_ * dh_;
    // Match PmlPartition.cpp: Simulation::c0_ / thickness_ * log10(1 / R_)
    // R is hardcoded to 1e-5 in PmlPartition.h/cpp
    double R_val = 1.0e-5;
    double zeta_max = Simulation::c0_ / thickness * log10(1.0 / R_val);
    std::vector<double> h_zetax(size, 0.0), h_zetay(size, 0.0), h_zetaz(size, 0.0);

    for (int k = 0; k < depth_; k++) {
        for (int j = 0; j < height_; j++) {
            for (int i = 0; i < width_; i++) {
                int idx = k * height_ * width_ + j * width_ + i;
                switch (type_) {
                case P_BACK:   h_zetaz[idx] = zeta_max * ((k + 1) * dh_ / thickness - sin(2 * M_PI * (k + 1) * dh_ / thickness) / 2 / M_PI); break;
                case P_FRONT:  h_zetaz[idx] = zeta_max * ((depth_ - k) * dh_ / thickness - sin(2 * M_PI * (depth_ - k) * dh_ / thickness) / 2 / M_PI); break;
                case P_BOTTOM: h_zetay[idx] = zeta_max * ((j + 1) * dh_ / thickness - sin(2 * M_PI * (j + 1) * dh_ / thickness) / 2 / M_PI); break;
                case P_TOP:    h_zetay[idx] = zeta_max * ((height_ - j) * dh_ / thickness - sin(2 * M_PI * (height_ - j) * dh_ / thickness) / 2 / M_PI); break;
                case P_RIGHT:  h_zetax[idx] = zeta_max * ((i + 1) * dh_ / thickness - sin(2 * M_PI * (i + 1) * dh_ / thickness) / 2 / M_PI); break;
                case P_LEFT:   h_zetax[idx] = zeta_max * ((width_ - i) * dh_ / thickness - sin(2 * M_PI * (width_ - i) * dh_ / thickness) / 2 / M_PI); break;
                }
            }
        }
    }
    CHECK_CUDA(cudaMemcpy(d_zetax_, h_zetax.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_zetay_, h_zetay.data(), bytes, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_zetaz_, h_zetaz.data(), bytes, cudaMemcpyHostToDevice));
}

PmlPartition::~PmlPartition()
{
    cudaFree(d_psi_); cudaFree(d_phi_x_); cudaFree(d_phi_y_); cudaFree(d_phi_z_);
    cudaFree(d_zetax_); cudaFree(d_zetay_); cudaFree(d_zetaz_);
    cudaFree(d_p_new_); cudaFree(d_v_new_); cudaFree(d_psi_new_);
    cudaFree(d_phi_x_new_); cudaFree(d_phi_y_new_); cudaFree(d_phi_z_new_);
}

void PmlPartition::Update()
{
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(
        (width_ + blockSize.x - 1) / blockSize.x, 
        (height_ + blockSize.y - 1) / blockSize.y, 
        (depth_ + blockSize.z - 1) / blockSize.z);

    PmlUpdateKernel<<<gridSize, blockSize>>>(
        d_pressure_, d_p_new_, d_velocity_, d_v_new_, 
        d_psi_, d_psi_new_,
        d_phi_x_, d_phi_x_new_,
        d_phi_y_, d_phi_y_new_,
        d_phi_z_, d_phi_z_new_,
        d_zetax_, d_zetay_, d_zetaz_,
        d_force_,
        width_, height_, depth_,
        dh_, Simulation::dt_, Simulation::c0_);
    
    // Swap buffers
    std::swap(d_pressure_, d_p_new_);
    std::swap(d_velocity_, d_v_new_);
    std::swap(d_psi_, d_psi_new_);
    std::swap(d_phi_x_, d_phi_x_new_);
    std::swap(d_phi_y_, d_phi_y_new_);
    std::swap(d_phi_z_, d_phi_z_new_);
}
