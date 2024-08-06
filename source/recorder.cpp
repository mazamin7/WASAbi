#include "recorder.h"
#include "simulation.h"
#include <iostream>


Recorder::Recorder(int x, int y, int z, int total_steps)
	: x_(x), y_(y), z_(z), total_steps_(total_steps)
{
	static int id_generator = 0;
	id_ = id_generator++;
	std::string filename;
	std::string dir_name = std::to_string(Simulation::dh_) + "_" + std::to_string(Partition::boundary_absorption_) + "_" + std::to_string(Simulation::air_absorption_alpha1_) + "_" + std::to_string(Simulation::air_absorption_alpha2_);
	filename = "./output/" + dir_name + "/out_" + std::to_string(id_) + ".txt";
	output_.open(filename, std::ios::out);
	filename = "./output/" + dir_name + "/response_" + std::to_string(id_) + ".txt";
	response_.open(filename, std::ios::out);
}


Recorder::~Recorder()
{
	output_.close();
	response_.close();
}

void Recorder::FindPartition(std::vector<std::shared_ptr<Partition>> partitions)
{
	for (auto partition : partitions)
	{
		if (partition->x_start_<x_ - 5 && partition->x_end_>x_ + 4 &&
			partition->y_start_<y_ - 5 && partition->y_end_>y_ + 4 &&
			partition->z_start_<z_ - 5 && partition->z_end_>z_ + 4)
		{
			part_ = partition;
			break;
		}
	}

	partitions_ = partitions;
}

void Recorder::RecordField(int time_step)
{
	bool render_pml = false;

	int x_start_, x_end_;
	int y_start_, y_end_;
	int z_start_, z_end_;

	int size_x_, size_y_, size_z_;

	x_start_ = y_start_ = z_start_ = std::numeric_limits<int>::max();
	x_end_ = y_end_ = z_end_ = std::numeric_limits<int>::min();

	for (auto partition : partitions_)
	{
		x_start_ = std::min(x_start_, partition->x_start_);
		y_start_ = std::min(y_start_, partition->y_start_);
		z_start_ = std::min(z_start_, partition->z_start_);

		x_end_ = std::max(x_end_, partition->x_end_);
		y_end_ = std::max(y_end_, partition->y_end_);
		z_end_ = std::max(z_end_, partition->z_end_);
	}

	size_x_ = x_end_ - x_start_;
	size_y_ = y_end_ - y_start_;
	size_z_ = z_end_ - z_start_;

	if ((time_step < total_steps_) && (time_step % 10 == 0))
	{
		double* values_ = (double*)calloc(size_x_ * size_y_ * size_z_, sizeof(double));

		for (auto partition : partitions_)
		{
			if (!render_pml)
			{
				if (!partition->should_render_) continue;
			}

			for (int i = 0; i < partition->width_; i++) {
				for (int j = 0; j < partition->height_; j++) {
					for (int k = 0; k < partition->depth_; k++) {
						values_[(partition->z_start_ + k) * size_y_ * size_x_ + (partition->y_start_ + j) * size_x_ + (partition->x_start_ + i)] = partition->get_pressure(i, j, k);
					}
				}
			}
		}

		for (int i = 0; i < size_x_ * size_y_ * size_z_; ++i) {
			output_ << values_[i] << " "; // Write each value followed by a newline
		}

		output_ << std::endl;
	}
}

void Recorder::RecordResponse(int time_step)
{
	if (time_step <= total_steps_)
	{
		response_ << part_->get_pressure(x_, y_, z_) << std::endl;
	}
}

std::vector<std::shared_ptr<Recorder>> Recorder::ImportRecorders(std::string path)
{
	std::vector<std::shared_ptr<Recorder>> recorders;

	std::ifstream file;
	file.open(path, std::ifstream::in);
	while (file.good())
	{
		int x, y, z;
		file >> x >> y >> z;
		if (file.eof()) break;
		recorders.push_back(std::make_shared<Recorder>(x / Simulation::dh_, y / Simulation::dh_, z / Simulation::dh_, Simulation::duration_ / Simulation::dt_));
	}
	file.close();

	return recorders;
}
