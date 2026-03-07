#include "cu_partition.h"
#include <cuda_runtime.h>
#include <stdexcept>
#include <stdio.h>
#include "sound_source.h"

#define CHECK_CUDA(call) \
    { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            printf("CUDA Error in " #call ": %s\n", cudaGetErrorString(err)); \
        } \
    }

CuPartition::CuPartition(int xs, int ys, int zs, int w, int h, int d)
    : Partition(xs, ys, zs, w, h, d)
{
    // Allocate core simulation arrays natively on the GPU VRAM
    size_t vol_size = width_ * height_ * depth_ * sizeof(double);
    cudaMalloc((void**)&d_pressure_, vol_size);
    cudaMalloc((void**)&d_velocity_, vol_size);
    cudaMalloc((void**)&d_force_, vol_size);
    cudaMalloc((void**)&d_residue_, vol_size);

    CHECK_CUDA(cudaMemset((void*)d_pressure_, 0, vol_size));
    CHECK_CUDA(cudaMemset((void*)d_velocity_, 0, vol_size));
    CHECK_CUDA(cudaMemset((void*)d_force_, 0, vol_size));
    CHECK_CUDA(cudaMemset((void*)d_residue_, 0, vol_size));

    // Initialize source buffers (assume max 100 sources per partition for now)
    num_sources_ = 0;
    cudaMalloc((void**)&d_source_indices_, 100 * sizeof(int));
    cudaMalloc((void**)&d_source_values_, 100 * sizeof(double));

    cudaDeviceSynchronize();
}

CuPartition::~CuPartition()
{
    cudaFree(d_pressure_);
    cudaFree(d_velocity_);
    cudaFree(d_force_);
    cudaFree(d_residue_);
    cudaFree(d_source_indices_);
    cudaFree(d_source_values_);
}

// ---- Helper: index into flat device array ----
static inline int idx3(int x, int y, int z, int w, int h)
{
    return z * h * w + y * w + x;
}

static double device_read(const double* d_arr, int offset)
{
    double val = 0.0;
    cudaError_t err = cudaMemcpy(&val, d_arr + offset, sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        printf("ERROR in device_read (offset %d): %s\n", offset, cudaGetErrorString(err));
    }
    return val;
}

static void device_write(double* d_arr, int offset, double val)
{
    cudaMemcpy(d_arr + offset, &val, sizeof(double), cudaMemcpyHostToDevice);
}

static void device_add(double* d_arr, int offset, double delta)
{
    double cur;
    cudaMemcpy(&cur, d_arr + offset, sizeof(double), cudaMemcpyDeviceToHost);
    cur += delta;
    cudaMemcpy(d_arr + offset, &cur, sizeof(double), cudaMemcpyHostToDevice);
}

// ---- pressure_field (returns null since data is on GPU) ----
double* CuPartition::get_pressure_field()
{
    return nullptr; // caller must use get_pressure() instead
}

// ---- Pressure accessors ----
double CuPartition::get_pressure(int x, int y, int z)
{
    return device_read(d_pressure_, idx3(x, y, z, width_, height_));
}
void CuPartition::set_pressure(int x, int y, int z, double v)
{
    device_write(d_pressure_, idx3(x, y, z, width_, height_), v);
}
void CuPartition::add_to_pressure(int x, int y, int z, double v)
{
    device_add(d_pressure_, idx3(x, y, z, width_, height_), v);
}

// ---- Velocity accessors ----
double CuPartition::get_velocity(int x, int y, int z)
{
    return device_read(d_velocity_, idx3(x, y, z, width_, height_));
}
void CuPartition::set_velocity(int x, int y, int z, double v)
{
    device_write(d_velocity_, idx3(x, y, z, width_, height_), v);
}
void CuPartition::add_to_velocity(int x, int y, int z, double v)
{
    device_add(d_velocity_, idx3(x, y, z, width_, height_), v);
}

// ---- Residue accessors ----
double CuPartition::get_residue(int x, int y, int z)
{
    return device_read(d_residue_, idx3(x, y, z, width_, height_));
}
void CuPartition::set_residue(int x, int y, int z, double v)
{
    device_write(d_residue_, idx3(x, y, z, width_, height_), v);
}
void CuPartition::add_to_residue(int x, int y, int z, double v)
{
    device_add(d_residue_, idx3(x, y, z, width_, height_), v);
}

// ---- Force accessors ----
double CuPartition::get_force(int x, int y, int z)
{
    return device_read(d_force_, idx3(x, y, z, width_, height_));
}
void CuPartition::set_force(int x, int y, int z, double v)
{
    device_write(d_force_, idx3(x, y, z, width_, height_), v);
}

// ---- Reset helpers ----
void CuPartition::reset_forces()
{
    cudaMemset((void*)d_force_, 0, width_ * height_ * depth_ * sizeof(double));
}

void CuPartition::reset_residues()
{
    cudaMemset((void*)d_residue_, 0, width_ * height_ * depth_ * sizeof(double));
}

void CuPartition::AddCuBoundary(CuBoundary* b)
{
    cu_boundaries_.push_back(b);
}

// ---- GPU Performance Overrides ----

__global__ void PostMergeKernel(
    double* d_velocity, const double* d_residue,
    int width, int height, int depth,
    double dt, double alpha1)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < width && y < height && z < depth) {
        int idx = (z * height * width) + (y * width) + x;
        double res = d_residue[idx];
        double update = (dt / (1.0 + 2.0 * dt * alpha1)) * res;
        d_velocity[idx] += update;
    }
}

void CuPartition::PostMerge()
{
    dim3 blockSize(8, 8, 8);
    dim3 gridSize(
        (width_ + blockSize.x - 1) / blockSize.x,
        (height_ + blockSize.y - 1) / blockSize.y,
        (depth_ + blockSize.z - 1) / blockSize.z);

    PostMergeKernel<<<gridSize, blockSize>>>(
        d_velocity_, d_residue_, width_, height_, depth_, dt_, air_absorption_alpha1_);
}

__global__ void ApplySourceKernel(double* d_force, const int* d_indices, const double* d_values, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        d_force[d_indices[i]] = d_values[i];
    }
}

void CuPartition::ComputeSourceForcingTerms(double t)
{
    if (sources_.empty()) return;

    if (num_sources_ != sources_.size()) {
        num_sources_ = (int)sources_.size();
        std::vector<int> h_indices;
        for (auto s : sources_) {
            h_indices.push_back(idx3(s->x() - x_start_, s->y() - y_start_, s->z() - z_start_, width_, height_));
        }
        cudaMemcpy(d_source_indices_, h_indices.data(), num_sources_ * sizeof(int), cudaMemcpyHostToDevice);
    }

    std::vector<double> h_values;
    for (auto s : sources_) {
        h_values.push_back(s->SampleValue(t));
    }
    cudaMemcpy(d_source_values_, h_values.data(), num_sources_ * sizeof(double), cudaMemcpyHostToDevice);

    int threads = 64;
    int blocks = (num_sources_ + threads - 1) / threads;
    ApplySourceKernel<<<blocks, threads>>>(d_force_, d_source_indices_, d_source_values_, num_sources_);
}

__device__ uint32_t RGBAToUint32(int r, int g, int b, int a) {
    return (uint32_t)((a << 24) | (b << 16) | (g << 8) | r);
}

__global__ void ColorMapKernel(
    const double* d_pressure, uint32_t* d_pixels,
    int width, int height, int depth,
    int plane_type, int coord,
    int screen_width, int screen_height,
    int x_offset, int y_offset, float v_coef, bool should_render)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    int p_w = 0, p_h = 0;
    if (plane_type == 0) { p_w = width; p_h = height; } // XY
    else if (plane_type == 1) { p_w = depth; p_h = height; } // YZ (looking from X)
    else if (plane_type == 2) { p_w = width; p_h = depth; } // XZ

    if (i < p_w && j < p_h) {
        int px = 0, py = 0, pz = 0;
        if (plane_type == 0) { px = i; py = j; pz = coord; }
        else if (plane_type == 1) { px = coord; py = j; pz = i; }
        else if (plane_type == 2) { px = i; py = coord; pz = j; }

        int idx = (pz * height * width) + (py * width) + px;
        double pressure = d_pressure[idx];
        double norm = 0.5 * fmax(-1.0, fmin(1.0, pressure * v_coef)) + 0.5;
        
        int r, g, b;
        if (norm >= 0.5) {
            r = (int)(255 - round(255.0 * 2.0 * (norm - 0.5)));
            g = (int)(255 - round(255.0 * 2.0 * (norm - 0.5)));
            b = 255;
        } else {
            r = 255;
            g = (int)(255 - round(255.0 * (1.0 - 2.0 * norm)));
            b = (int)(255 - round(255.0 * (1.0 - 2.0 * norm)));
        }

        if (!should_render) {
            r *= 0.5; g *= 0.5; b *= 0.5;
        }

        int out_x = x_offset + i;
        int out_y = y_offset + j;
        if (out_x >= 0 && out_x < screen_width && out_y >= 0 && out_y < screen_height) {
            d_pixels[out_y * screen_width + out_x] = RGBAToUint32(r, g, b, 255);
        }
    }
}

void CuPartition::RenderToBuffer(uint32_t* d_pixels, int plane_type, int coord, int screen_width, int screen_height, int x_offset, int y_offset, float v_coef)
{
    int p_w = 0, p_h = 0, local_coord = 0;
    if (plane_type == 0) { p_w = width_; p_h = height_; local_coord = coord - z_start_; if (local_coord < 0 || local_coord >= depth_) return; }
    else if (plane_type == 1) { p_w = depth_; p_h = height_; local_coord = coord - x_start_; if (local_coord < 0 || local_coord >= width_) return; }
    else if (plane_type == 2) { p_w = width_; p_h = depth_; local_coord = coord - y_start_; if (local_coord < 0 || local_coord >= height_) return; }

    dim3 blockSize(16, 16);
    dim3 gridSize((p_w + blockSize.x - 1) / blockSize.x, (p_h + blockSize.y - 1) / blockSize.y);

    ColorMapKernel<<<gridSize, blockSize>>>(
        d_pressure_, d_pixels, width_, height_, depth_,
        plane_type, local_coord, screen_width, screen_height,
        x_offset, y_offset, v_coef, should_render_);
}

std::vector<double> CuPartition::get_xy_plane(int z)
{
    int local_z = z - z_start_;
    int slice_size = width_ * height_;
    std::vector<double> plane(slice_size);
    
    if (local_z >= 0 && local_z < depth_) {
        size_t offset = local_z * slice_size;
        cudaMemcpy(plane.data(), d_pressure_ + offset, slice_size * sizeof(double), cudaMemcpyDeviceToHost);
    }
    return plane;
}

std::vector<double> CuPartition::get_yz_plane(int x)
{
    // YZ plane is non-contiguous in our layout (Z, Y, X)
    // We could write a kernel to gather it, but for now we fallback 
    // to the base class if it's not the primary visualization axis.
    // However, if we want high speed, we should gather it.
    return Partition::get_yz_plane(x); 
}

std::vector<double> CuPartition::get_xz_plane(int y)
{
    return Partition::get_xz_plane(y);
}
