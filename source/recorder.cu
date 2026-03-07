#include "recorder.h"
#include "simulation.h"
#include <iostream>
#include <fstream>
#include <string>


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
	
	// NOTE THAT IT SAVES 1 EVERY 10 TIME STEPS
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
    
    // 1. Attempt to open the file
    file.open(path, std::ifstream::in);

    // --- CHECK 1: File Open Failure ---
    if (!file.is_open()) 
    {
        std::cerr << "FATAL ERROR: Failed to open recorder import file." << std::endl;
        std::cerr << "  Path attempted: " << path << std::endl;
        // Since we can't proceed without the input, we return an empty list.
        return recorders; 
    }
    std::cout << "SUCCESS: Recorder import file opened: " << path << std::endl;

    int line_count = 0;
    while (true)
    {
        int x, y, z;
        
        // 2. Attempt to read x, y, z values
        file >> x >> y >> z;
        line_count++;

        // --- CHECK 2: Data Read Failure (Bad Format) ---
        // If the stream is not in a 'good' state after the read attempt (e.g., hit non-numeric text),
        // but it's not the end-of-file yet, we report an error and break.
        if (file.fail() && !file.eof()) 
        {
            std::cerr << "ERROR: Invalid data format detected in recorder file." << std::endl;
            std::cerr << "  Problem occurred near line: " << line_count << std::endl;
            break; 
        }

        // Check 3: End of file was reached *after* attempting to read, 
        // or a read failure occurred right at the end of the file.
        if (file.eof()) break; 
        
        // --- CHECK 4: Sanity Check (Preventing Division by Zero) ---
        // Assuming Simulation::dh_ and Simulation::dt_ must be non-zero for calculations
        if (Simulation::dh_ == 0.0 || Simulation::dt_ == 0.0)
        {
             std::cerr << "FATAL ERROR: Simulation parameters (dh_ or dt_) are zero. Cannot calculate grid coordinates." << std::endl;
             file.close();
             return std::vector<std::shared_ptr<Recorder>>();
        }
        
        // If all checks pass, create the Recorder object
        recorders.push_back(std::make_shared<Recorder>((int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_), (int)(Simulation::duration_ / Simulation::dt_)));
    }
    
    file.close();

    if (recorders.empty() && line_count > 0)
    {
        std::cerr << "WARNING: File was read, but no valid recorders were created (check Simulation parameters)." << std::endl;
    }
    else if (!recorders.empty())
    {
        std::cout << "SUCCESS: Imported " << recorders.size() << " recorders." << std::endl;
    }

    return recorders;
}
