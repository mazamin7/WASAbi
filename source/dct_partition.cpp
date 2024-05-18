#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"


DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d)
	: Partition(xs, ys, zs, w, h, d), pressure_(w, h, d), velocity_(w, h, d), residue_(w, h, d), force_(w, h, d)
{
	should_render_ = true;
	info_.type = "DCT";

	next_pressure_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	next_velocity_modes_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

	w0_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

	alpha_ = air_absorption_;

	lx2_ = width_ * width_ * dh_ * dh_;
	ly2_ = height_ * height_ * dh_ * dh_;
	lz2_ = depth_ * depth_ * dh_ * dh_;

	for (int i = 1; i <= depth_; i++)
	{
		for (int j = 1; j <= height_; j++)
		{
			for (int k = 1; k <= width_; k++)
			{
				int idx = (i - 1) * height_ * width_ + (j - 1) * width_ + (k - 1);
				double w0 = c0_ * M_PI * sqrt((i - 1) * (i - 1) / lz2_ + (j - 1) * (j - 1) / ly2_ + (k - 1) * (k - 1) / lx2_);
				w0_[idx] = w0;
			}
		}
	}
}


DctPartition::~DctPartition()
{
	free(next_pressure_modes_);
	free(next_velocity_modes_);
	free(w0_);
}

void DctPartition::Update()
{
	pressure_.ExecuteDct(); // current pressure
	velocity_.ExecuteDct(); // current pressure velocity
	force_.ExecuteDct(); // current force

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				int idx = i * height_ * width_ + j * width_ + k;

				if ((idx == 0) && (alpha_ == 0.0))
				{
					next_velocity_modes_[idx] = velocity_.modes_[idx] + dt_ * force_.modes_[idx];
					next_pressure_modes_[idx] = dt_ * velocity_.modes_[idx] + pressure_.modes_[idx] + dt_*dt_/2 * force_.modes_[idx];
				}
				else if ((idx == 0) && (alpha_ > 0.0))
				{
					double e2at = exp(-2 * alpha_ * dt_);

					next_velocity_modes_[idx] = e2at * velocity_.modes_[idx] + (1 - e2at)/(2 * alpha_) * force_.modes_[idx];
					next_pressure_modes_[idx] = (1 - e2at)/(2 * alpha_) * velocity_.modes_[idx] + pressure_.modes_[idx] + (e2at - 1) / (4 * alpha_ * alpha_) * force_.modes_[idx];
				}
				else if ((idx > 0) && (alpha_ < w0_[idx]))
				{
					double inv_w02 = 1 / w0_[idx] / w0_[idx];
					double alpha_sqr = alpha_ * alpha_;

					double w = sqrt(w0_[idx] * w0_[idx] - alpha_sqr);

					double cwt = cos(w * dt_);
					double swt = sin(w * dt_);
					double eatm = exp(-alpha_ * dt_);

					double w2 = w * w;
					double inv_w = 1 / w;

					double xe = force_.modes_[idx] * inv_w02;

					next_velocity_modes_[idx] = eatm * (velocity_.modes_[idx] * (cwt - alpha_ * inv_w * swt) - (w + alpha_sqr * inv_w) * (pressure_.modes_[idx] - xe) * swt);
					next_pressure_modes_[idx] = xe + eatm * ((pressure_.modes_[idx] - xe) * (cwt + alpha_ * inv_w * swt) + swt * inv_w * velocity_.modes_[idx]);
				}
				else if ((idx > 0) && (alpha_ > w0_[idx]))
				{
					double inv_w02 = 1 / w0_[idx] / w0_[idx];
					double alpha_sqr = alpha_ * alpha_;

					double alphad = sqrt(alpha_sqr - w0_[idx] * w0_[idx]);
					double alpha1 = alpha_ + alphad;
					double alpha2 = alpha_ - alphad;

					double eat1 = exp(-alpha1 * dt_);
					double eat2 = exp(-alpha2 * dt_);

					next_velocity_modes_[idx] = (0.5*(eat1 + eat2) + 0.5/alphad * alpha_ * (eat1 - eat2)) * velocity_.modes_[idx] + (-0.5 * (alpha1 * eat1 +  alpha2 * eat2) - 0.5 / alphad * alpha_ * (alpha2 * eat2 - alpha1 * eat1)) * pressure_.modes_[idx] + inv_w02 * 0.5 * (alpha1 * eat1 + alpha2 * eat2 + alpha_/alphad * (alpha2 * eat2 - alpha1 * eat1)) * force_.modes_[idx];
					next_pressure_modes_[idx] = 0.5 / alphad * (eat2 - eat1) * velocity_.modes_[idx] + (eat1 + eat2 + 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1)) * pressure_.modes_[idx] + inv_w02 * (1 - alpha1 - alpha2 - 0.5 / alphad * (alpha2 * eat2 - alpha1 * eat1)) * force_.modes_[idx];
				}
			}
		}
	}

	memcpy((void*)velocity_.modes_, (void*)next_velocity_modes_, depth_ * width_ * height_ * sizeof(double));
	memcpy((void*)pressure_.modes_, (void*)next_pressure_modes_, depth_ * width_ * height_ * sizeof(double));

	velocity_.ExecuteIdct();
	pressure_.ExecuteIdct();
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

void DctPartition::reset_forces()
{
	force_.reset();
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
