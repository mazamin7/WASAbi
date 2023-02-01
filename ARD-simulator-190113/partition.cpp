#include "partition.h"
#include "boundary.h"
#include "sound_source.h"
#include "dct_partition.h"
#include "simulation.h"
#include "tools.h"
#include <fstream>
#include <iostream>

/*
The constructor creates a new Partition object, initializing its member variables.

The member variables x_start_, y_start_, z_start_, width_, height_, and depth_ are 
initialized with the values passed in as parameters xs, ys, zs, w, h, and d, respectively.

Next, id_generator is a static variable that is incremented each time a new Partition 
object is created. The current value of id_generator is assigned to info_.id.

The member variables dh_, dt_, and c0_ are assigned the values of the static variables 
Simulation::dh_, Simulation::dt_, and Simulation::c0_, respectively.

The end coordinates of the partition along each axis (x_end_, y_end_, and z_end_) are 
calculated by adding the width, height, and depth of the partition to the start 
coordinates along each axis.

Finally, three nested loops are used to initialize the member variables 
left_free_borders_, right_free_borders_, top_free_borders_, bottom_free_borders_, 
front_free_borders_, and back_free_borders_ to vectors of vectors of integers, each 
with the value true. 

The dimensions of each vector depend on the depth, height, and width of the partition.
*/
Partition::Partition(int xs, int ys, int zs, int w, int h, int d)
	: x_start_(xs), y_start_(ys), z_start_(zs), width_(w), height_(h), depth_(d)
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

/*
AddBoundary() adds a Boundary object to a Partition object. 

It updates the number of boundaries associated with the Partition object, and 
updates the boolean values in the arrays left_free_borders_, 
right_free_borders_, top_free_borders_, bottom_free_borders_, indicating 
whether there is a boundary at each location. 

The function uses a Boundary object's type 
(Boundary::X_BOUNDARY or Boundary::Y_BOUNDARY) and start/end positions 
(x_start_, x_end_, y_start_, y_end_) to determine the location of the boundary 
and update the corresponding left_free_borders_, right_free_borders_, 
top_free_borders_, or bottom_free_borders_ array.
*/
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

/*
The function ImportPartitions is a member function of the class Partition. 

It reads partition information from a file located at the path path and returns 
a vector of shared pointers to Partition objects, which represent the partitions.

The function starts by opening a file stream file to the specified path. 

It then repeatedly reads the start coordinates (x_start, y_start, z_start) and 
the dimensions (width, height, depth) of a partition until the end of the file 
is reached. 

For each partition, the function creates a new shared pointer to a DctPartition 
object, initializing it with the start coordinates and dimensions scaled by the 
Simulation class constant dh_. 

The created partition is then added to the end of the partitions vector. 

Finally, the file stream is closed and the partitions vector is returned.
*/
std::vector<std::shared_ptr<Partition>> Partition::ImportPartitions(std::string path)
{
	std::vector<std::shared_ptr<Partition>> partitions;

	std::ifstream file;
	file.open(path, std::ifstream::in);
	while (file.good())
	{
		int x_start, y_start, z_start;
		int width, height, depth;
		file >> x_start >> y_start >> z_start;
		file >> width >> height >> depth;
		if (file.eof()) break;

		partitions.push_back(std::make_shared<DctPartition>(x_start / Simulation::dh_, y_start / Simulation::dh_, z_start / Simulation::dh_, width / Simulation::dh_, height / Simulation::dh_, depth / Simulation::dh_));
	}
	file.close();
	return partitions;
}

void Partition::Info()
{
	std::cout << info_.type << " Partition "<< info_.id <<": "
		<< x_start_ << "," << y_start_ << "," << z_start_ << "->"
		<< x_end_ << "," << y_end_ << "," << z_end_ << std::endl;
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
