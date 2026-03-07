#include "boundary.h"
#include "partition.h"
#include "simulation.h"
#include <algorithm>


Boundary::Boundary(BoundaryType type, double absorp, std::shared_ptr<Partition> a, std::shared_ptr<Partition> b,
	int xs, int xe, int ys, int ye, int zs, int ze)
	: type_(type), boundary_absorption_(absorp), a_(a), b_(b), x_start_(xs), x_end_(xe), y_start_(ys), y_end_(ye), z_start_(zs), z_end_(ze)
{
	static int id_generator = 0;
	info_.id = id_generator++;
	info_.a_id = a_->info_.id;
	info_.b_id = b_->info_.id;
}

Boundary::~Boundary()
{
}

__constant__ double d_coefs[6][6] = {
    {  0.0,   0.0,   -2.0,    2.0,   0.0,  0.0 },
    {  0.0,  -2.0,   27.0,  -27.0,   2.0,  0.0 },
    { -2.0,  27.0, -270.0,  270.0, -27.0,  2.0 },
    {  2.0, -27.0,  270.0, -270.0,  27.0, -2.0 },
    {  0.0,   2.0,  -27.0,   27.0,  -2.0,  0.0 },
    {  0.0,   0.0,    2.0,   -2.0,   0.0,  0.0 } 
};

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

void Boundary::ComputeResidues()
{
    dim3 blockSize(16, 16);
    auto a_cu = a_;
    auto b_cu = b_;

    if (type_ == X_BOUNDARY)
    {
        bool is_a_left = (x_start_ <= a_cu->x_end_ && x_end_ >= a_cu->x_end_);
        auto left = is_a_left ? a_cu : b_cu;
        auto right = is_a_left ? b_cu : a_cu;

        int y_len = y_end_ - y_start_;
        int z_len = z_end_ - z_start_;
        dim3 gridSize((z_len + blockSize.x - 1) / blockSize.x, (y_len + blockSize.y - 1) / blockSize.y);

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

//void Boundary::ComputeForcingTerms()	// Micheal Oliver edition.
//{
//	double coefs[][7] = {
//					{ 0,   0,    0,   0,    0,   -2,  2 },
//					{ 0,   0,    0,   -2,   27,  -27, 2 },
//					{ 0,   -2,   27,  -270, 270, -27, 2 },
//					{ 2,   -27,  270, -270, 27,  -2,  0 },
//					{ 2,   -27,  27,  -2,   0,   0,   0 },
//					{ 2,   -2,   0,   0,    0,   0,   0 } };
//	if (type_ == X_BOUNDARY)
//	{
//		bool is_a_left = (x_start_ <= a_->x_end_&&x_end_ >= a_->x_end_);
//		auto left = is_a_left ? a_ : b_;
//		auto right = is_a_left ? b_ : a_;
//		for (int i = y_start_; i < y_end_; i++)
//		{
//			for (int j = z_start_; j < z_end_; j++)
//			{
//				int left_y = i - left->y_start_;
//				int right_y = i - right->y_start_;
//
//				int left_z = j - left->z_start_;
//				int right_z = j - right->z_start_;
//				
//				for (int m = -3; m < 3; m++)
//				{
//					int left_x = left->width_ + m;
//					int right_x = m;
//					double sip = 0.0;
//					double sip1 = 0.0;
//					double sip2 = 0.0;
//					double fi = 0.0;
//
//					int coefs_idx = 0;
//
//					for (int k = left_x - 3; k < left->width_; k++, coefs_idx++)
//					{
//						sip1 = coefs[m+3][coefs_idx] * left->get_pressure(k, left_y, left_z);
//					}
//
//					for (int k = 0; coefs_idx < 7; k++, coefs_idx++)
//					{
//						sip2 = coefs[m+3][coefs_idx] * right->get_pressure(k, right_y, right_z);
//					}
//
//					if (m < 0)
//					{
//						sip = sip1 + sip2;
//						fi = sip * (Simulation::c0_*Simulation::c0_) / (180.0*Simulation::dh_*Simulation::dh_);
//						left->set_force(left_x, left_y, left_z, fi);
//					}
//					else {
//						sip = sip1 + sip2;
//						fi = sip * (Simulation::c0_*Simulation::c0_) / (180.0*Simulation::dh_*Simulation::dh_);
//						right->set_force(right_x, right_y, right_z, fi);
//					}
//				}
//
//			}
//		}
//	}
//	else if (type_ == Y_BOUNDARY)
//	{
//
//	}
//}


std::shared_ptr<Boundary> Boundary::FindBoundary(std::shared_ptr<Partition> a, std::shared_ptr<Partition> b, double absorp)
{
	int xa_min = a->x_start_;
	int xa_max = xa_min + a->width_;
	int xb_min = b->x_start_;
	int xb_max = xb_min + b->width_;
	int x_overlapped = std::min(xa_max, xb_max) - std::max(xa_min, xb_min);

	int ya_min = a->y_start_;
	int ya_max = ya_min + a->height_;
	int yb_min = b->y_start_;
	int yb_max = yb_min + b->height_;
	int y_overlapped = std::min(ya_max, yb_max) - std::max(ya_min, yb_min);

	int za_min = a->z_start_;
	int za_max = za_min + a->depth_;
	int zb_min = b->z_start_;
	int zb_max = zb_min + b->depth_;
	int z_overlapped = std::min(za_max, zb_max) - std::max(za_min, zb_min);

	bool both_cuda = true; // All partitions are now CUDA partitions

	if (x_overlapped == 0 && y_overlapped > 0 && z_overlapped > 0)
	{
		bool is_right_boundary = (xa_max == xb_min);
		int x_start = (is_right_boundary ? xa_max - 3 : xb_max - 3);
		int x_end = x_start + 6;
		int y_start = std::max(ya_min, yb_min);
		int y_end = y_start + y_overlapped;
		int z_start = std::max(za_min, zb_min);
		int z_end = z_start + z_overlapped;

		if (both_cuda)
			return std::make_shared<Boundary>(X_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
		else
			return std::make_shared<Boundary>(X_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
	}
	else if (y_overlapped == 0 && x_overlapped > 0 && z_overlapped > 0)
	{
		bool is_bottom_boundary = (ya_max == yb_min);
		int x_start = std::max(xa_min, xb_min);
		int x_end = x_start + x_overlapped;
		int y_start = (is_bottom_boundary ? ya_max - 3 : yb_max - 3);
		int y_end = y_start + 6;
		int z_start = std::max(za_min, zb_min);
		int z_end = z_start + z_overlapped;

		if (both_cuda)
			return std::make_shared<Boundary>(Y_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
		else
			return std::make_shared<Boundary>(Y_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
	}
	else if (z_overlapped == 0 && x_overlapped > 0 && y_overlapped > 0)
	{
		bool is_back_boundary = (za_max == zb_min);
		int x_start = std::max(xa_min, xb_min);
		int x_end = x_start + x_overlapped;
		int y_start = std::max(ya_min, yb_min);
		int y_end = y_start + y_overlapped;
		int z_start = (is_back_boundary ? za_max - 3 : zb_max - 3);
		int z_end = z_start + 6;

		if (both_cuda)
			return std::make_shared<Boundary>(Z_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
		else
			return std::make_shared<Boundary>(Z_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
	}
	return nullptr;
}

void Boundary::Info()
{
	switch (type_)
	{
	case X_BOUNDARY:
		std::cout << "Boundary " << info_.id << ": X " << info_.a_id << "," << info_.b_id << " | " << x_start_ << "," << y_start_ << "," << z_start_ << "->" << x_end_ << "," << y_end_ << "," << z_end_ << std::endl;
		break;
	case Y_BOUNDARY:
		std::cout << "Boundary " << info_.id << ": Y " << info_.a_id << "," << info_.b_id << " | " << x_start_ << "," << y_start_ << "," << z_start_ << "->" << x_end_ << "," << y_end_ << "," << z_end_ << std::endl;
		break;
	case Z_BOUNDARY:
		std::cout << "Boundary " << info_.id << ": Z " << info_.a_id << "," << info_.b_id << " | " << x_start_ << "," << y_start_ << "," << z_start_ << "->" << x_end_ << "," << y_end_ << "," << z_end_ << std::endl;
		break;
	default:
		break;
	}
}
