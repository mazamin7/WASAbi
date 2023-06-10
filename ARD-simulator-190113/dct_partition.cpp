#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"


DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d, double alpha_abs)
	: Partition(xs, ys, zs, w, h, d, alpha_abs), pressure_(w, h, d), velocity_(w, h, d), residue_(w, h, d), force_(w, h, d), force_r_(w, h, d)
{
	should_render_ = true;
	info_.type = "DCT";

	second_order_ = alpha_abs == 0;

	prev_pressure_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	next_pressure_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	prev_velocity_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	next_velocity_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

	cwt_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	swt_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	w_omega_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	w2_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	inv_w_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	inv_w2_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

	alpha_ = alpha_abs;
	alpha2_ = alpha_ * alpha_;
	eatm_ = exp(-alpha_ * dt_);

	lx2_ = width_ * width_*dh_*dh_;
	ly2_ = height_ * height_*dh_*dh_;
	lz2_ = depth_ * depth_*dh_*dh_;

	for (int i = 1; i <= depth_; i++)
	{
		for (int j = 1; j <= height_; j++)
		{
			for (int k = 1; k <= width_; k++)
			{
				int idx = (i - 1) * height_ * width_ + (j - 1) * width_ + (k - 1);
				double w = c0_ * M_PI * sqrt((i - 1) * (i - 1) / lz2_ + (j - 1) * (j - 1) / ly2_ + (k - 1) * (k - 1) / lx2_);
				cwt_[idx] = cos(w * dt_);
				swt_[idx] = sin(w * dt_);
				w_omega_[idx] = w;
				w2_[idx] = w * w;
				inv_w_[idx] = 1 / w;
				inv_w2_[idx] = inv_w_[idx] * inv_w_[idx];
			}
		}
	}
}


DctPartition::~DctPartition()
{
	free(prev_pressure_modes_);
	free(next_pressure_modes_);
	free(prev_velocity_modes_);
	free(next_velocity_modes_);
	free(cwt_);
	free(swt_);
	free(w_omega_);
	free(w2_);
	free(inv_w_);
	free(inv_w2_);
}

void DctPartition::Update_pressure()
{
	// prev_pressure_modes_ previous pressure
	pressure_.ExecuteDct(); // current pressure
	velocity_.ExecuteDct(); // current velocity
	force_.ExecuteDct(); // current force
	residue_.ExecuteDct(); // current residual

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				int idx = i * height_ * width_ + j * width_ + k;
				
				if (idx == 0) {
					if (second_order_)
						// p_next_dct(n) = 2 * p_curr_dct(n) - p_prev_dct(n) + dt*dt * force_dct(n);
						next_pressure_modes_[idx] = (1 - 1e-10) * ( 2.0 * pressure_.modes_[idx] - prev_pressure_modes_[idx] + dt_ * dt_ * (force_.modes_[idx] + residue_.modes_[idx]));
					else{
						// p_next_dct(n) = p_curr_dct(n) + dt * v_curr_dct(n);
						next_pressure_modes_[idx] = (1 - 1e-10) * ( pressure_.modes_[idx] + dt_ * velocity_.modes_[idx] );
					}
				}else{
					if (second_order_)
						// p_next_dct(n) = 2 * p_curr_dct(n) .* cwt(n) - p_prev_dct(n) + (2 * force_dct(n) . / w2(n)) .* (1 - cwt(n));
						next_pressure_modes_[idx] = (1 - 1e-10) * ( 2.0 * pressure_.modes_[idx] * cwt_[idx] - prev_pressure_modes_[idx] + (2.0 * (force_.modes_[idx] + residue_.modes_[idx]) * inv_w2_[idx]) * (1.0 - cwt_[idx]) );
					else{
						double xe = force_.modes_[idx] * inv_w2_[idx];
						// p_next_dct(n) = xe + eatm * ((p_curr_dct(n) - xe) .* (cwt(n) + alpha_abs * inv_w(n) .* swt(n)) + swt(n) .* inv_w(n) .* v_curr_dct(n));
						next_pressure_modes_[idx] = (1 - 1e-10) * (xe + eatm_ * ((pressure_.modes_[idx] - xe) * (cwt_[idx] + alpha_ * inv_w_[idx] * swt_[idx]) + swt_[idx] * inv_w_[idx] * velocity_.modes_[idx]) );
					}
				}
			}
		}
	}

	memcpy((void *)prev_pressure_modes_, (void *)pressure_.modes_, depth_ * width_ * height_ * sizeof(double));
	memcpy((void *)pressure_.modes_, (void *)next_pressure_modes_, depth_ * width_ * height_ * sizeof(double));

	pressure_.ExecuteIdct();
}

void DctPartition::Update_velocity()
{
	if (second_order_) // velocity is not considered in the second order case
		return;

	// prev_pressure_modes_ current pressure
	velocity_.ExecuteDct(); // current velocity
	force_.ExecuteDct(); // next force
	force_r_.ExecuteDct(); // next corrected force

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				int idx = i * height_ * width_ + j * width_ + k;

				if (idx == 0) {
					// v_next_dct(n) = (v_curr_dct(n) + dt * force_dct(n)) / (1 + 2 * dt * alpha_abs);
					next_velocity_modes_[idx] = (velocity_.modes_[idx] + dt_ * force_r_.modes_[idx]) / (1 + 2 * dt_ * alpha_abs_);
				}
				else {
					double xe = force_r_.modes_[idx] * inv_w2_[idx];
					// v_next_dct(n) = eatm * (v_curr_dct(n) .* (cwt(n) - alpha_abs * inv_w(n) .* swt(n)) - (w(n) + alpha2 * inv_w(n)) .* (p_curr_dct(n) - xe) .* swt(n));
					next_velocity_modes_[idx] = eatm_ * (velocity_.modes_[idx] * (cwt_[idx] - alpha_ * inv_w_[idx] * swt_[idx]) - (w_omega_[idx] + alpha2_ * inv_w_[idx]) * (prev_pressure_modes_[idx] - xe) * swt_[idx]);
				}
			}
		}
	}

	memcpy((void*)prev_velocity_modes_, (void*)velocity_.modes_, depth_ * width_ * height_ * sizeof(double));
	memcpy((void*)velocity_.modes_, (void*)next_velocity_modes_, depth_ * width_ * height_ * sizeof(double));

	velocity_.ExecuteIdct();
}

double* DctPartition::get_pressure_field()
{
	return pressure_.values_;
}

double DctPartition::get_pressure(int x, int y, int z)
{
	return pressure_.get_value(x, y, z);
}

void DctPartition::set_pressure(int x, int y, int z, double v)
{
	pressure_.set_value(x, y, z, v);
}

void DctPartition::add_to_pressure(int x, int y, int z, double v)
{
	pressure_.set_value(x, y, z, pressure_.get_value(x, y, z) + v);
}

double DctPartition::get_velocity(int x, int y, int z)
{
	return velocity_.get_value(x, y, z);
}

void DctPartition::set_velocity(int x, int y, int z, double v)
{
	velocity_.set_value(x, y, z, v);
}

void DctPartition::add_to_velocity(int x, int y, int z, double v)
{
	velocity_.set_value(x, y, z, velocity_.get_value(x, y, z) + v);
}

double DctPartition::get_residue(int x, int y, int z)
{
	return residue_.get_value(x, y, z);
}

void DctPartition::set_residue(int x, int y, int z, double v)
{
	residue_.set_value(x, y, z, v);
}

void DctPartition::add_to_residue(int x, int y, int z, double v)
{
	residue_.set_value(x, y, z, residue_.get_value(x, y, z) + v);
}

double DctPartition::get_force(int x, int y, int z)
{
	return force_.get_value(x, y, z);
}

void DctPartition::set_force(int x, int y, int z, double v)
{
	force_.set_value(x, y, z, v);
}

double DctPartition::get_force_r(int x, int y, int z)
{
	return force_r_.get_value(x, y, z);
}

void DctPartition::set_force_r(int x, int y, int z, double v)
{
	force_r_.set_value(x, y, z, v);
}

void DctPartition::reset_forces()
{
	force_.reset();
	force_r_.reset();
}

void DctPartition::reset_residues()
{
	residue_.reset();
}

std::vector<double> DctPartition::get_xy_forcing_plane(int z)
{
	std::vector<double> xy_plane;
	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			xy_plane.push_back(get_force(j, i, z));
		}
	}
	return xy_plane;
}

std::vector<double> DctPartition::get_xy_force_plane(int z)
{
	std::vector<double> xy_plane;
	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			xy_plane.push_back(get_force(j, i, z));
		}
	}
	return xy_plane;
}
