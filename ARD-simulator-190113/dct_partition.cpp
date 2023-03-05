#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"


DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d, double h_abs)
	: Partition(xs, ys, zs, w, h, d, h_abs), pressure_(w, h, d), force_(w, h, d), residue_(w, h, d)
{
	should_render_ = true;
	info_.type = "DCT";

	is_damped_ = h_abs != 0;

	prev_modes_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	next_modes_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	velocity_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

	cwt_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	swt_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	alpha_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	alpha2_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	eatm_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	w_omega_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	inv_w_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));
	inv_w2_ = (double*)calloc(width_ * height_ * depth_, sizeof(double));

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
				double w = c0_ * M_PI * sqrt(i * i / lz2_ + j * j / ly2_ + k * k / lx2_);
				cwt_[idx] = cos(w * dt_);
				swt_[idx] = sin(w * dt_);
				alpha_[idx] = w * h_abs;
				alpha2_[idx] = alpha_[idx] * alpha_[idx];
				eatm_[idx] = exp(-alpha_[idx] * dt_);
				w_omega_[idx] = w;
				inv_w_[idx] = 1 / w;
				inv_w2_[idx] = inv_w_[idx] * inv_w_[idx];
			}
		}
	}
}


DctPartition::~DctPartition()
{
	free(prev_modes_);
	free(next_modes_);
	free(velocity_);
	free(cwt_);
	free(swt_);
	free(alpha_);
	free(alpha2_);
	free(eatm_);
	free(w_omega_);
	free(inv_w_);
	free(inv_w2_);
}

void DctPartition::Update()
{
	force_.ExcuteDct();
	pressure_.ExcuteDct();

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				int idx = i * height_ * width_ + j * width_ + k;
				
				if (!is_damped_)
					next_modes_[idx] = 2.0 * pressure_.modes_[idx] * cwt_[idx] - prev_modes_[idx] + (2.0 * force_.modes_[idx] * inv_w2_[idx]) * (1.0 - cwt_[idx]);
				else{
					double xe = force_.modes_[idx] * inv_w2_[idx];
					next_modes_[idx] = xe + eatm_[idx] * ((pressure_.modes_[idx] - xe) * (cwt_[idx] + alpha_[idx] * inv_w_[idx] * swt_[idx]) + swt_[idx] * inv_w_[idx] * velocity_[idx]);
					velocity_[idx] = eatm_[idx] * (velocity_[idx] * (cwt_[idx] - alpha_[idx] * inv_w_[idx] * swt_[idx]) - (w_omega_[idx] + alpha2_[idx] * inv_w_[idx]) * (pressure_.modes_[idx] - xe) * swt_[idx]);
				}
			}
		}
	}

	memcpy((void *)prev_modes_, (void *)pressure_.modes_, depth_ * width_ * height_ * sizeof(double));
	memcpy((void *)pressure_.modes_, (void *)next_modes_, depth_ * width_ * height_ * sizeof(double));

	pressure_.ExcuteIdct();
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

double DctPartition::get_residue(int x, int y, int z)
{
	return residue_.get_value(x, y, z);
}

void DctPartition::set_residue(int x, int y, int z, double v)
{
	residue_.set_value(x, y, z, v);
}

double DctPartition::get_force(int x, int y, int z)
{
	return force_.get_value(x, y, z);
}

void DctPartition::set_force(int x, int y, int z, double f)
{
	force_.set_value(x, y, z, f);
}

void DctPartition::reset_forces()
{
	force_.reset();
}

void DctPartition::reset_residues()
{
	residue_.reset();
}

double DctPartition::check_reset_residues()
{
	double acc = 0.0;

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				double temp = get_residue(k, j, i);
				acc = acc + temp;
			}
		}
	}

	return acc;
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
