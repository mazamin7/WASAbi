#include "partition.h"
#include "boundary.h"
#include "sound_source.h"
#include "dct_partition.h"
#include "simulation.h"
#include "tools.h"
#include <fstream>
#include <iostream>

Partition::Partition(int xs, int ys, int zs, int w, int h, int d, double h_abs)
	: x_start_(xs), y_start_(ys), z_start_(zs), width_(w), height_(h), depth_(d), h_abs_(h_abs)
{
	static int id_generator = 0;
	info_.id = id_generator++;
	dh_ = Simulation::dh_;
	dt_ = Simulation::dt_;
	c0_ = Simulation::c0_;
	x_end_ = x_start_ + width_;
	y_end_ = y_start_ + height_;
	z_end_ = z_start_ + depth_;

	for (int i = 0; i < depth_; i++)
	{
		std::vector<int> v;
		for (int j = 0; j < height_; j++)
		{
			v.push_back(true);
		}
		left_free_borders_.push_back(v);
		right_free_borders_.push_back(v);
	}
	for (int i = 0; i < depth_; i++)
	{
		std::vector<int> v;
		for (int j = 0; j < width_; j++)
		{
			v.push_back(true);
		}
		top_free_borders_.push_back(v);
		bottom_free_borders_.push_back(v);
	}
	for (int i = 0; i < height_; i++)
	{
		std::vector<int> v;
		for (int j = 0; j < width_; j++)
		{
			v.push_back(true);
		}
		front_free_borders_.push_back(v);
		back_free_borders_.push_back(v);
	}
}

Partition::~Partition()
{
}

std::vector<double> Partition::get_xy_plane(int z)
{
	std::vector<double> xy_plane;
	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			xy_plane.push_back(get_pressure(j, i, z));
		}
	}
	return xy_plane;
}

std::vector<double> Partition::get_yz_plane(int x)
{
	std::vector<double> yz_plane;
	for (int i = 0; i < depth_; i++) {
		for (int j = 0; j < height_; j++) {
			yz_plane.push_back(get_pressure(x, j, i));
		}
	}
	return yz_plane;
}

std::vector<double> Partition::get_xz_plane(int y)
{
	std::vector<double> xz_plane;
	for (int i = 0; i < depth_; i++) {
		for (int j = 0; j < width_; j++) {
			xz_plane.push_back(get_pressure(j, y, i));
		}
	}
	return xz_plane;
}

std::vector<double> Partition::get_xy_forcing_plane(int z)
{
	return std::vector<double>();
}

void Partition::AddBoundary(std::shared_ptr<Boundary> boundary)
{
	info_.num_boundaries++;
	if (boundary->type_ == Boundary::X_BOUNDARY)
	{
		if (boundary->x_start_ < x_start_ && boundary->x_end_ > x_start_)	// left boundary
		{
			for (int i = 0; i < depth_; i++)
			{
				for (int j = boundary->y_start_; j < boundary->y_end_; j++)
				{
					left_free_borders_[i][j - y_start_] = false;
				}
			}
		}
		else {
			for (int i = 0; i < depth_; i++)
			{
				for (int j = boundary->y_start_; j < boundary->y_end_; j++)
				{
					right_free_borders_[i][j - y_start_] = false;
				}
			}
		}
	}
	else if (boundary->type_ == Boundary::Y_BOUNDARY)
	{
		if (boundary->y_start_ < y_start_ && boundary->y_end_ > y_start_)		// top boundary
		{
			for (int i = 0; i < depth_; i++)
			{
				for (int j = boundary->x_start_; j < boundary->x_end_; j++)
				{
					top_free_borders_[i][j - x_start_] = false;
				}
			}
		}
		else {
			for (int i = 0; i < depth_; i++)
			{
				for (int j = boundary->x_start_; j < boundary->x_end_; j++)
				{
					bottom_free_borders_[i][j - x_start_] = false;
				}
			}
		}
	}
}

void Partition::AddSource(std::shared_ptr<SoundSource> source)
{
	info_.num_sources++;
	sources_.push_back(source);
}

std::vector<std::shared_ptr<Partition>> Partition::ImportPartitions(std::string path)
{
	std::vector<std::shared_ptr<Partition>> partitions;

	std::ifstream file;
	file.open(path, std::ifstream::in);
	while (file.good())
	{
		int x_start, y_start, z_start;
		int width, height, depth;
		double h_abs;

		file >> x_start >> y_start >> z_start;
		file >> width >> height >> depth;
		file >> h_abs;

		if (file.eof()) break;

		partitions.push_back(std::make_shared<DctPartition>(x_start / Simulation::dh_, y_start / Simulation::dh_, z_start / Simulation::dh_, width / Simulation::dh_, height / Simulation::dh_, depth / Simulation::dh_, h_abs));
	}
	file.close();
	return partitions;
}

void Partition::Info()
{
	std::cout << info_.type << " Partition " << info_.id << ": "
		<< x_start_ << "," << y_start_ << "," << z_start_ << "->"
		<< x_end_ << "," << y_end_ << "," << z_end_ << "; "
		<< "Dimensionless abs. coeff.: " << h_abs_ << std::endl;
	std::cout << "    -> " << std::to_string(info_.num_sources) << " sources; "
		<< std::to_string(info_.num_boundaries) << " boundaries; " << std::endl;
}

void Partition::ComputeSourceForcingTerms(double t)
{
	for (auto source : sources_)
	{
		set_force(
			source->x_ - x_start_,
			source->y_ - y_start_,
			source->z_ - z_start_,
			source->SampleValue(t));
	}
}

void Partition::PreMerge()
{
	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				auto old = get_force(k, j, i);
				auto res = get_residue(k, j, i);
				set_force(k, j, i, res + old);
			}
		}
	}
}

void Partition::PostMerge()
{
	double dt2 = Simulation::dt_ * Simulation::dt_;

	for (int i = 0; i < depth_; i++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int k = 0; k < width_; k++)
			{
				auto old = get_pressure(k, j, i);
				auto res = get_residue(k, j, i);
				set_pressure(k, j, i, res * dt2 + old);
			}
		}
	}
}
