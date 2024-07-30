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

void Boundary::ComputeResidues()
{
	double coefs[6][6] = {
		{  0.0,   0.0,   -2.0,    2.0,   0.0,  0.0 },
		{  0.0,  -2.0,   27.0,  -27.0,   2.0,  0.0 },
		{ -2.0,  27.0, -270.0,  270.0, -27.0,  2.0 },
		{  2.0, -27.0,  270.0, -270.0,  27.0, -2.0 },
		{  0.0,   2.0,  -27.0,   27.0,  -2.0,  0.0 },
		{  0.0,   0.0,    2.0,   -2.0,   0.0,  0.0 } };

	if (type_ == X_BOUNDARY)
	{
		bool is_a_left = (x_start_ <= a_->x_end_ && x_end_ >= a_->x_end_);
		auto left = is_a_left ? a_ : b_;
		auto right = is_a_left ? b_ : a_;
//#pragma omp parallel for
		for (int i = y_start_; i < y_end_; i++)
		{
			for (int j = z_start_; j < z_end_; j++)
			{
				int left_y = i - left->y_start_;
				int right_y = i - right->y_start_;
				int left_z = j - left->z_start_;
				int right_z = j - right->z_start_;
				int left_x, right_x;
				double sip, sip1, sip2, res;

				// Unrolled iteration for m = -3
				left_x = left->width_ - 3;
				right_x = -3;

				sip1 = coefs[0][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[0][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[0][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[0][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[0][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[0][5] * right->get_pressure(2, right_y, right_z);

				sip = left->include_self_terms_ ? (sip1 + sip2) : sip2;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				left->add_to_residue(left_x, left_y, left_z, boundary_absorption_ * res);

				// Unrolled iteration for m = -2
				left_x = left->width_ - 2;
				right_x = -2;

				sip1 = coefs[1][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[1][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[1][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[1][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[1][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[1][5] * right->get_pressure(2, right_y, right_z);

				sip = left->include_self_terms_ ? (sip1 + sip2) : sip2;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				left->add_to_residue(left_x, left_y, left_z, boundary_absorption_ * res);

				// Unrolled iteration for m = -1
				left_x = left->width_ - 1;
				right_x = -1;

				sip1 = coefs[2][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[2][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[2][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[2][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[2][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[2][5] * right->get_pressure(2, right_y, right_z);

				sip = left->include_self_terms_ ? (sip1 + sip2) : sip2;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				left->add_to_residue(left_x, left_y, left_z, boundary_absorption_ * res);

				// Unrolled iteration for m = 0
				left_x = left->width_;
				right_x = 0;

				sip1 = coefs[3][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[3][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[3][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[3][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[3][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[3][5] * right->get_pressure(2, right_y, right_z);

				sip = right->include_self_terms_ ? (sip1 + sip2) : sip1;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				right->add_to_residue(right_x, right_y, right_z, boundary_absorption_ * res);

				// Unrolled iteration for m = 1
				left_x = left->width_ + 1;
				right_x = 1;

				sip1 = coefs[4][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[4][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[4][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[4][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[4][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[4][5] * right->get_pressure(2, right_y, right_z);

				sip = right->include_self_terms_ ? (sip1 + sip2) : sip1;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				right->add_to_residue(right_x, right_y, right_z, boundary_absorption_ * res);

				// Unrolled iteration for m = 2
				left_x = left->width_ + 2;
				right_x = 2;

				sip1 = coefs[5][0] * left->get_pressure(left->width_ - 3, left_y, left_z);
				sip1 += coefs[5][1] * left->get_pressure(left->width_ - 2, left_y, left_z);
				sip1 += coefs[5][2] * left->get_pressure(left->width_ - 1, left_y, left_z);

				sip2 = coefs[5][3] * right->get_pressure(0, right_y, right_z);
				sip2 += coefs[5][4] * right->get_pressure(1, right_y, right_z);
				sip2 += coefs[5][5] * right->get_pressure(2, right_y, right_z);

				sip = right->include_self_terms_ ? (sip1 + sip2) : sip1;
				res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
				right->add_to_residue(right_x, right_y, right_z, boundary_absorption_ * res);

			}
		}
	}
	else if (type_ == Y_BOUNDARY)
	{
		bool is_a_top = (y_start_ <= a_->y_end_ && y_end_ >= a_->y_end_);
		auto top = is_a_top ? a_ : b_;
		auto bottom = is_a_top ? b_ : a_;
// #pragma omp parallel for
		for (int i = x_start_; i < x_end_; i++)
		{
			for (int j = z_start_; j < z_end_; j++)
			{
				int top_x = i - top->x_start_;
				int bottom_x = i - bottom->x_start_;
				int top_z = j - top->z_start_;
				int bottom_z = j - bottom->z_start_;
				// Unrolled iteration for m = -3
				{
					int top_y = top->height_ - 3;
					int bottom_y = -3;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[0][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[0][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[0][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[0][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[0][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[0][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = top->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					top->add_to_residue(top_x, top_y, top_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = -2
				{
					int top_y = top->height_ - 2;
					int bottom_y = -2;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[1][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[1][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[1][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[1][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[1][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[1][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = top->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					top->add_to_residue(top_x, top_y, top_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = -1
				{
					int top_y = top->height_ - 1;
					int bottom_y = -1;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[2][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[2][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[2][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[2][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[2][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[2][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = top->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					top->add_to_residue(top_x, top_y, top_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 0
				{
					int top_y = top->height_;
					int bottom_y = 0;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[3][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[3][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[3][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[3][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[3][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[3][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = bottom->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					bottom->add_to_residue(bottom_x, bottom_y, bottom_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 1
				{
					int top_y = top->height_ + 1;
					int bottom_y = 1;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[4][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[4][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[4][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[4][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[4][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[4][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = bottom->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					bottom->add_to_residue(bottom_x, bottom_y, bottom_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 2
				{
					int top_y = top->height_ + 2;
					int bottom_y = 2;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[5][0] * top->get_pressure(top_x, top->height_ - 3 + 0, top_z);
					sip1 += coefs[5][1] * top->get_pressure(top_x, top->height_ - 3 + 1, top_z);
					sip1 += coefs[5][2] * top->get_pressure(top_x, top->height_ - 3 + 2, top_z);

					sip2 += coefs[5][3] * bottom->get_pressure(bottom_x, 0, bottom_z);
					sip2 += coefs[5][4] * bottom->get_pressure(bottom_x, 1, bottom_z);
					sip2 += coefs[5][5] * bottom->get_pressure(bottom_x, 2, bottom_z);

					sip = bottom->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					bottom->add_to_residue(bottom_x, bottom_y, bottom_z, boundary_absorption_ * res);
				}

			}
		}
	}
	else if (type_ == Z_BOUNDARY)
	{
		bool is_a_front = (z_start_ <= a_->z_end_ && z_end_ >= a_->z_end_);
		auto front = is_a_front ? a_ : b_;
		auto back = is_a_front ? b_ : a_;
// #pragma omp parallel for
		for (int i = x_start_; i < x_end_; i++)
		{
			for (int j = y_start_; j < y_end_; j++)
			{
				int front_x = i - front->x_start_;
				int back_x = i - back->x_start_;
				int front_y = j - front->y_start_;
				int back_y = j - back->y_start_;
				// Unrolled iteration for m = -3
				{
					int front_z = front->depth_ - 3;
					int back_z = -3;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[0][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[0][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[0][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[0][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[0][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[0][5] * back->get_pressure(back_x, back_y, 2);

					sip = front->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					front->add_to_residue(front_x, front_y, front_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = -2
				{
					int front_z = front->depth_ - 2;
					int back_z = -2;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[1][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[1][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[1][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[1][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[1][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[1][5] * back->get_pressure(back_x, back_y, 2);

					sip = front->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					front->add_to_residue(front_x, front_y, front_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = -1
				{
					int front_z = front->depth_ - 1;
					int back_z = -1;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[2][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[2][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[2][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[2][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[2][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[2][5] * back->get_pressure(back_x, back_y, 2);

					sip = front->include_self_terms_ ? (sip1 + sip2) : sip2;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					front->add_to_residue(front_x, front_y, front_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 0
				{
					int front_z = front->depth_;
					int back_z = 0;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[3][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[3][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[3][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[3][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[3][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[3][5] * back->get_pressure(back_x, back_y, 2);

					sip = back->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					back->add_to_residue(back_x, back_y, back_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 1
				{
					int front_z = front->depth_ + 1;
					int back_z = 1;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[4][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[4][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[4][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[4][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[4][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[4][5] * back->get_pressure(back_x, back_y, 2);

					sip = back->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					back->add_to_residue(back_x, back_y, back_z, boundary_absorption_ * res);
				}

				// Unrolled iteration for m = 2
				{
					int front_z = front->depth_ + 2;
					int back_z = 2;

					double sip1 = 0.0;
					double sip2 = 0.0;
					double sip;
					double res;

					sip1 += coefs[5][0] * front->get_pressure(front_x, front_y, front->depth_ - 3);
					sip1 += coefs[5][1] * front->get_pressure(front_x, front_y, front->depth_ - 2);
					sip1 += coefs[5][2] * front->get_pressure(front_x, front_y, front->depth_ - 1);

					sip2 += coefs[5][3] * back->get_pressure(back_x, back_y, 0);
					sip2 += coefs[5][4] * back->get_pressure(back_x, back_y, 1);
					sip2 += coefs[5][5] * back->get_pressure(back_x, back_y, 2);

					sip = back->include_self_terms_ ? (sip1 + sip2) : sip1;
					res = sip * (Simulation::c0_ * Simulation::c0_) / (180.0 * Simulation::dh_ * Simulation::dh_);
					back->add_to_residue(back_x, back_y, back_z, boundary_absorption_ * res);
				}

			}
		}
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
	// > 0 ==> partitions overlapping
	// = 0 ==> partitions share exact boundary
	// < 0 ==> partitions not touching

	int ya_min = a->y_start_;
	int ya_max = ya_min + a->height_;
	int yb_min = b->y_start_;
	int yb_max = yb_min + b->height_;
	int y_overlapped = std::min(ya_max, yb_max) - std::max(ya_min, yb_min);

	int z_start = a->z_start_;
	int z_end = z_start + a->depth_;

	if (x_overlapped == 0 && y_overlapped > 0)
	{
		bool is_right_boundary = (xa_max == xb_min);

		int x_start = (is_right_boundary ? xa_max - 3 : xb_max - 3);
		int x_end = x_start + 6;

		int y_start = ya_min + std::max(0, yb_min - ya_min);
		int y_end = y_start + y_overlapped;

		return std::make_shared<Boundary>(X_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
	}
	else if (y_overlapped == 0 && x_overlapped > 0)
	{
		bool is_bottom_boundary = (ya_max == yb_min);

		int x_start = xa_min + std::max(0, xb_min - xa_min);
		int x_end = x_start + x_overlapped;

		int y_start = (is_bottom_boundary ? ya_max - 3 : yb_max - 3);
		int y_end = y_start + 6;

		return std::make_shared<Boundary>(Y_BOUNDARY, absorp, a, b, x_start, x_end, y_start, y_end, z_start, z_end);
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
