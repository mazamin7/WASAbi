#define _USE_MATH_DEFINES
#include <cmath>
#include "dct_partition.h"

/*
The constructor sets the values of various data members including pressure_, force_, 
prev_modes_, next_modes_, cwt_, w2_, and others.

It first initializes the base class Partition with the given arguments 
xs, ys, zs, w, h, and d which represents the starting position xs, ys, zs and the 
width w, height h, and depth d of the partition.

Then it sets the value of should_render_ to true and sets the type of the info_ member 
to "DCT".

The member prev_modes_, next_modes_, cwt_, and w2_ are dynamically allocated memory of 
size width_ * height_ * depth_ using the calloc function.

The member lx2_, ly2_, and lz2_ are computed based on the values of width_, height_, 
and depth_.

Finally, a nested for loop is used to compute values of w2_ (squared angular frequency) 
and cwt_ ( cos(omega*dt) ).

The loop iterates over the entire volume of the partition defined by width_, height_, 
and depth_ and computes the values of w2_ and cwt_ at each point. 

The computation is based on the values of lx2_, ly2_, and lz2_ and the constant 
c0_ (speed of sound) and dt_ (time step).
*/
DctPartition::DctPartition(int xs, int ys, int zs, int w, int h, int d)
	: Partition(xs, ys, zs, w, h, d), pressure_(w, h, d), force_(w, h, d)
{
	should_render_ = true;
	info_.type = "DCT";

	prev_modes_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	next_modes_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	cwt_ = (double*)calloc(width_*height_*depth_, sizeof(double));
	w2_ = (double*)calloc(width_*height_*depth_, sizeof(double));

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
				w2_[idx] = w * w;
				cwt_[idx] = cos(w * dt_);
			}
		}
	}
}


DctPartition::~DctPartition()
{
	free(prev_modes_);
	free(next_modes_);
	free(cwt_);
	free(w2_);
}

/*
The Update function of the DctPartition class performs the following steps:

1. The force_ object executes the discrete cosine transform (DCT) with the 
ExcuteDct() function.

2. For each cell in the partition, the value of next_modes_[idx] is updated 
based on the current values of pressure_.modes_[idx], prev_modes_[idx], and 
force_.modes_[idx].

3. The value of prev_modes_ is updated with the current values of 
pressure_.modes_.

4. The values of pressure_.modes_ are updated with the values of next_modes_.

5. The pressure_ object executes the inverse discrete cosine transform (IDCT) 
with the ExcuteIdct function.
*/
void DctPartition::Update()
{
	force_.ExcuteDct();
	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				int idx = i * height_ * width_ + j * width_ + k;
				next_modes_[idx] = 0.999 * (2.0 * pressure_.modes_[idx] * cwt_[idx] - prev_modes_[idx] + (2.0 * force_.modes_[idx] / w2_[idx]) * (1.0 - cwt_[idx]));
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

void DctPartition::set_force(int x, int y, int z, double f)
{
	force_.set_value(x, y, z, f);
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

double DctPartition::get_force(int x, int y, int z)
{
	return force_.get_value(x, y, z);
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
