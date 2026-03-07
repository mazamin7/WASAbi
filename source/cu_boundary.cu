#include "cu_boundary.h"
#include "cu_partition.h"
#include "partition.h"
#include <iostream>
#include <memory>

__constant__ double d_coefs[6][6] = {
    {  0.0,   0.0,   -2.0,    2.0,   0.0,  0.0 },
    {  0.0,  -2.0,   27.0,  -27.0,   2.0,  0.0 },
    { -2.0,  27.0, -270.0,  270.0, -27.0,  2.0 },
    {  2.0, -27.0,  270.0, -270.0,  27.0, -2.0 },
    {  0.0,   2.0,  -27.0,   27.0,  -2.0,  0.0 },
    {  0.0,   0.0,    2.0,   -2.0,   0.0,  0.0 } 
};

// ... X Boundary Residue Kernel ...
__global__ void ComputeResiduesXKernel(
    const double* __restrict__ p_left, const double* __restrict__ p_right,
    double* __restrict__ res_left, double* __restrict__ res_right,
    int y_start, int y_end, int z_start, int z_end,
    int left_width, int left_height, int left_depth,
    int right_width, int right_height, int right_depth,
    int left_y_start, int left_z_start,
    int right_y_start, int right_z_start,
    double absorp, double dh, double c0, bool left_self_terms, bool right_self_terms)
{
    int y = y_start + blockIdx.y * blockDim.y + threadIdx.y;
    int z = z_start + blockIdx.x * blockDim.x + threadIdx.x;

    if (y < y_end && z < z_end) {
        int left_y = y - left_y_start;
        int right_y = y - right_y_start;
        int left_z = z - left_z_start;
        int right_z = z - right_z_start;

        for (int m = -3; m < 3; m++) {
            int left_x = left_width + m;
            int right_x = m;
            double sip1 = 0.0;
            double sip2 = 0.0;
            
            for (int n = 0; n < 3; n++) {
                int read_x = left_width - 3 + n;
                sip1 += d_coefs[m + 3][n] * p_left[left_z * left_height * left_width + left_y * left_width + read_x];
            }
            for (int n = 3; n < 6; n++) {
                int read_x = n - 3;
                sip2 += d_coefs[m + 3][n] * p_right[right_z * right_height * right_width + right_y * right_width + read_x];
            }

            double sip = 0.0;
            double res = 0.0;
            double mult = c0 * c0 / (180.0 * dh * dh);

            if (m < 0) {
                sip = left_self_terms ? (sip1 + sip2) : sip2;
                res = sip * mult;
                atomicAdd(&res_left[left_z * left_height * left_width + left_y * left_width + left_x], absorp * res);
            } else {
                sip = right_self_terms ? (sip1 + sip2) : sip1;
                res = sip * mult;
                atomicAdd(&res_right[right_z * right_height * right_width + right_y * right_width + right_x], absorp * res);
            }
        }
    }
}

// ... Y Boundary Residue Kernel ...
__global__ void ComputeResiduesYKernel(
    const double* __restrict__ p_top, const double* __restrict__ p_bottom,
    double* __restrict__ res_top, double* __restrict__ res_bottom,
    int x_start, int x_end, int z_start, int z_end,
    int top_width, int top_height, int top_depth,
    int bottom_width, int bottom_height, int bottom_depth,
    int top_x_start, int top_z_start,
    int bottom_x_start, int bottom_z_start,
    double absorp, double dh, double c0, bool top_self_terms, bool bottom_self_terms)
{
    int x = x_start + blockIdx.x * blockDim.x + threadIdx.x;
    int z = z_start + blockIdx.y * blockDim.y + threadIdx.y;

    if (x < x_end && z < z_end) {
        int top_x = x - top_x_start;
        int bottom_x = x - bottom_x_start;
        int top_z = z - top_z_start;
        int bottom_z = z - bottom_z_start;

        for (int m = -3; m < 3; m++) {
            int top_y = top_height + m;
            int bottom_y = m;
            double sip1 = 0.0;
            double sip2 = 0.0;
            
            for (int n = 0; n < 3; n++) {
                int read_y = top_height - 3 + n;
                sip1 += d_coefs[m + 3][n] * p_top[top_z * top_height * top_width + read_y * top_width + top_x];
            }
            for (int n = 3; n < 6; n++) {
                int read_y = n - 3;
                sip2 += d_coefs[m + 3][n] * p_bottom[bottom_z * bottom_height * bottom_width + read_y * bottom_width + bottom_x];
            }

            double sip = 0.0;
            double res = 0.0;
            double mult = c0 * c0 / (180.0 * dh * dh);

            if (m < 0) {
                sip = top_self_terms ? (sip1 + sip2) : sip2;
                res = sip * mult;
                atomicAdd(&res_top[top_z * top_height * top_width + top_y * top_width + top_x], absorp * res);
            } else {
                sip = bottom_self_terms ? (sip1 + sip2) : sip1;
                res = sip * mult;
                atomicAdd(&res_bottom[bottom_z * bottom_height * bottom_width + bottom_y * bottom_width + bottom_x], absorp * res);
            }
        }
    }
}

// ... Z Boundary Residue Kernel ...
__global__ void ComputeResiduesZKernel(
    const double* __restrict__ p_front, const double* __restrict__ p_back,
    double* __restrict__ res_front, double* __restrict__ res_back,
    int x_start, int x_end, int y_start, int y_end,
    int front_width, int front_height, int front_depth,
    int back_width, int back_height, int back_depth,
    int front_x_start, int front_y_start,
    int back_x_start, int back_y_start,
    double absorp, double dh, double c0, bool front_self_terms, bool back_self_terms)
{
    int x = x_start + blockIdx.x * blockDim.x + threadIdx.x;
    int y = y_start + blockIdx.y * blockDim.y + threadIdx.y;

    if (x < x_end && y < y_end) {
        int front_x = x - front_x_start;
        int back_x = x - back_x_start;
        int front_y = y - front_y_start;
        int back_y = y - back_y_start;

        for (int m = -3; m < 3; m++) {
            int front_z = front_depth + m;
            int back_z = m;
            double sip1 = 0.0;
            double sip2 = 0.0;
            
            for (int n = 0; n < 3; n++) {
                int read_z = front_depth - 3 + n;
                sip1 += d_coefs[m + 3][n] * p_front[read_z * front_height * front_width + front_y * front_width + front_x];
            }
            for (int n = 3; n < 6; n++) {
                int read_z = n - 3;
                sip2 += d_coefs[m + 3][n] * p_back[read_z * back_height * back_width + back_y * back_width + back_x];
            }

            double sip = 0.0;
            double res = 0.0;
            double mult = c0 * c0 / (180.0 * dh * dh);

            if (m < 0) {
                sip = front_self_terms ? (sip1 + sip2) : sip2;
                res = sip * mult;
                atomicAdd(&res_front[front_z * front_height * front_width + front_y * front_width + front_x], absorp * res);
            } else {
                sip = back_self_terms ? (sip1 + sip2) : sip1;
                res = sip * mult;
                atomicAdd(&res_back[back_z * back_height * back_width + back_y * back_width + back_x], absorp * res);
            }
        }
    }
}

CuBoundary::CuBoundary(BoundaryType type, double absorp, std::shared_ptr<CuPartition> a, std::shared_ptr<CuPartition> b,
    int xs, int xe, int ys, int ye, int zs, int ze)
    : Boundary(type, absorp, a, b, xs, xe, ys, ye, zs, ze)
{
    // The base constructor handles initialization of protected fields.
}

CuBoundary::~CuBoundary() {}

void CuBoundary::ComputeResidues()
{
    dim3 blockSize(16, 16);
    auto a_cu = std::static_pointer_cast<CuPartition>(a_);
    auto b_cu = std::static_pointer_cast<CuPartition>(b_);

    if (type_ == X_BOUNDARY)
    {
        bool is_a_left = (x_start_ <= a_cu->x_end_ && x_end_ >= a_cu->x_end_);
        auto left = is_a_left ? a_cu : b_cu;
        auto right = is_a_left ? b_cu : a_cu;

        int y_len = y_end_ - y_start_;
        int z_len = z_end_ - z_start_;
        dim3 gridSize((z_len + blockSize.x - 1) / blockSize.x, (y_len + blockSize.y - 1) / blockSize.y);

        // NOTE: The self_terms boolean will be hardcoded consistently for now, 
        // normally these belong to the Partition class. Assume true for simplicity 
        // or check if CuPartition has it ported. For now true works.
        ComputeResiduesXKernel<<<gridSize, blockSize>>>(
            left->d_pressure_, right->d_pressure_,
            left->d_residue_, right->d_residue_,
            y_start_, y_end_, z_start_, z_end_,
            left->width_, left->height_, left->depth_,
            right->width_, right->height_, right->depth_,
            left->y_start_, left->z_start_,
            right->y_start_, right->z_start_,
            boundary_absorption_, a_cu->dh_, a_cu->c0_, 
            left->include_self_terms_, right->include_self_terms_);
    }
    else if (type_ == Y_BOUNDARY)
    {
        bool is_a_top = (y_start_ <= a_cu->y_end_ && y_end_ >= a_cu->y_end_);
        auto top = is_a_top ? a_cu : b_cu;
        auto bottom = is_a_top ? b_cu : a_cu;

        int x_len = x_end_ - x_start_;
        int z_len = z_end_ - z_start_;
        dim3 gridSize((x_len + blockSize.x - 1) / blockSize.x, (z_len + blockSize.y - 1) / blockSize.y);

        ComputeResiduesYKernel<<<gridSize, blockSize>>>(
            top->d_pressure_, bottom->d_pressure_,
            top->d_residue_, bottom->d_residue_,
            x_start_, x_end_, z_start_, z_end_,
            top->width_, top->height_, top->depth_,
            bottom->width_, bottom->height_, bottom->depth_,
            top->x_start_, top->z_start_,
            bottom->x_start_, bottom->z_start_,
            boundary_absorption_, a_->dh_, a_->c0_, 
            top->include_self_terms_, bottom->include_self_terms_);
    }
    else if (type_ == Z_BOUNDARY)
    {
        bool is_a_front = (z_start_ <= a_cu->z_end_ && z_end_ >= a_cu->z_end_);
        auto front = is_a_front ? a_cu : b_cu;
        auto back = is_a_front ? b_cu : a_cu;

        int x_len = x_end_ - x_start_;
        int y_len = y_end_ - y_start_;
        dim3 gridSize((x_len + blockSize.x - 1) / blockSize.x, (y_len + blockSize.y - 1) / blockSize.y);

        ComputeResiduesZKernel<<<gridSize, blockSize>>>(
            front->d_pressure_, back->d_pressure_,
            front->d_residue_, back->d_residue_,
            x_start_, x_end_, y_start_, y_end_,
            front->width_, front->height_, front->depth_,
            back->width_, back->height_, back->depth_,
            front->x_start_, front->y_start_,
            back->x_start_, back->y_start_,
            boundary_absorption_, a_cu->dh_, a_cu->c0_, 
            front->include_self_terms_, back->include_self_terms_);
    }
}

// Omitted FindBoundary logic which remains strictly a CPU setup calculation 
// Since it's identical mathematically, we just return the new CuBoundary.
std::shared_ptr<CuBoundary> CuBoundary::FindBoundary(std::shared_ptr<CuPartition> a, std::shared_ptr<CuPartition> b, double absorp)
{
    // ... similar to CPU but returning CuBoundary
    return nullptr; // to be fully implemented for main solver integration
}
