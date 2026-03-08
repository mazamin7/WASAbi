#include "recorder.h"
#include "simulation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include "json.hpp"

using json = nlohmann::json;

Recorder::Recorder(const Config& config, int x, int y, int z, int total_steps, std::string dir_path)
	: config_(config), x_(x), y_(y), z_(z), total_steps_(total_steps), dir_path_(dir_path)
{
	static int id_generator = 0;
	id_ = id_generator++;
	
	output_path_ = dir_path_ + "/record_data_" + std::to_string(id_) + ".bin";
	response_path_ = dir_path_ + "/response_data_" + std::to_string(id_) + ".bin";

	cudaMalloc((void**)&d_response_buffer_, total_steps_ * sizeof(double));
    response_count_ = 0;
}

Recorder::~Recorder()
{
	if (field_buffer_) free(field_buffer_);
    if (d_response_buffer_) cudaFree(d_response_buffer_);
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

	// Calculate global simulation dimensions once
	go_x_ = go_y_ = go_z_ = std::numeric_limits<int>::max();
	int xe = std::numeric_limits<int>::min();
	int ye = std::numeric_limits<int>::min();
	int ze = std::numeric_limits<int>::min();

	for (auto p : partitions_) {
		go_x_ = std::min(go_x_, p->x_start_);
		go_y_ = std::min(go_y_, p->y_start_);
		go_z_ = std::min(go_z_, p->z_start_);
		xe = std::max(xe, p->x_end_);
		ye = std::max(ye, p->y_end_);
		ze = std::max(ze, p->z_end_);
	}
	gs_x_ = xe - go_x_;
	gs_y_ = ye - go_y_;
	gs_z_ = ze - go_z_;

	// Pre-allocate field buffer to avoid mid-loop overhead
	if (field_buffer_) free(field_buffer_);
	field_buffer_ = (double*)calloc(gs_x_ * gs_y_ * gs_z_, sizeof(double));
}

void Recorder::RecordField(int time_step)
{
	bool render_pml = false;

	// Record EVERY time step to match expected frame count
	if (time_step < total_steps_)
	{
		if (!output_.is_open()) {
			output_.open(output_path_, std::ios::out | std::ios::binary);
		}
		
		// Clear local host buffer before assembling
		memset(field_buffer_, 0, gs_x_ * gs_y_ * gs_z_ * sizeof(double));

		for (auto partition : partitions_)
		{
			if (!render_pml)
			{
				if (!partition->should_render_) continue;
			}

			// OPTIMIZATION: Bulk transfer from GPU to CPU per partition 
			// Instead of million individual get_pressure() calls
			size_t vol_size = partition->width_ * partition->height_ * partition->depth_ * sizeof(double);
			double* h_partition_field = (double*)malloc(vol_size);
			cudaMemcpy(h_partition_field, partition->d_pressure_, vol_size, cudaMemcpyDeviceToHost);

			// Assemble the partition field into the global field buffer
			for (int k = 0; k < partition->depth_; k++) {
				for (int j = 0; j < partition->height_; j++) {
					// Single memcpy for the entire row (X-dimension)
					size_t row_bytes = partition->width_ * sizeof(double);
					size_t src_offset = (k * partition->height_ * partition->width_) + (j * partition->width_);
					size_t dst_offset = ((partition->z_start_ - go_z_ + k) * gs_y_ * gs_x_) + 
					                   ((partition->y_start_ - go_y_ + j) * gs_x_) + 
					                   (partition->x_start_ - go_x_);
					
					memcpy(field_buffer_ + dst_offset, h_partition_field + src_offset, row_bytes);
				}
			}
			free(h_partition_field);
		}

		// Dump directly as binary
        int total_points = gs_x_ * gs_y_ * gs_z_;
		output_.write(reinterpret_cast<const char*>(field_buffer_), total_points * sizeof(double));
	}
}

__global__ void RecordResponseKernel(const double* d_pressure, double* d_buffer, int idx, int x, int y, int z, int w, int h)
{
    d_buffer[idx] = d_pressure[z * h * w + y * w + x];
}

void Recorder::RecordResponse(int time_step)
{
	if (time_step < total_steps_)
	{
		// FIX: Launch a tiny kernel to record value asynchronously on the GPU
        RecordResponseKernel<<<1, 1, 0, part_->stream_>>>(
            part_->d_pressure_, d_response_buffer_, response_count_++,
            x_ - part_->x_start_, y_ - part_->y_start_, z_ - part_->z_start_,
            part_->width_, part_->height_
        );
	}
}

void Recorder::FlushResponse()
{
	if (response_count_ == 0) return;

    std::vector<double> h_buffer(response_count_);
    cudaMemcpy(h_buffer.data(), d_response_buffer_, response_count_ * sizeof(double), cudaMemcpyDeviceToHost);

	if (!response_.is_open()) {
		response_.open(response_path_, std::ios::out | std::ios::binary);
	}
	response_.write(reinterpret_cast<const char*>(h_buffer.data()), response_count_ * sizeof(double));
	response_.close();
}

std::vector<std::shared_ptr<Recorder>> Recorder::ImportRecorders(const Config& config, std::string path, int total_steps__, std::string dir_path)
{
	std::vector<std::shared_ptr<Recorder>> recorders;

    if (std::filesystem::path(path).extension() == ".json") {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "WARNING: Could not open JSON asset file: " << path << std::endl;
            return recorders;
        }
        try {
            json j;
            file >> j;
            if (j.contains("recorders")) {
                for (auto& r : j["recorders"]) {
                    recorders.push_back(std::make_shared<Recorder>(config, 
                        (int)((double)r["x"] / Simulation::dh_), 
                        (int)((double)r["y"] / Simulation::dh_), 
                        (int)((double)r["z"] / Simulation::dh_), 
                        total_steps__, dir_path));
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "ERROR parsing JSON asset file: " << e.what() << std::endl;
        }
        return recorders;
    }

	std::ifstream file(path);
	if (!file.is_open())
	{
		std::cerr << "WARNING: Could not open unified asset file: " << path << std::endl;
		return recorders;
	}

	std::string line;
	while (std::getline(file, line))
	{
		if (line.empty()) continue;

		std::stringstream ss(line);
		std::string first_token;
		ss >> first_token;

		if (first_token == "R" || first_token == "r") {
			// It's explicitly a recorder
			double x, y, z, dummy;
			if ((ss >> x >> y >> z) && !(ss >> dummy)) {
				recorders.push_back(std::make_shared<Recorder>(config, (int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_), total_steps__, dir_path));
			}
		}
		else {
			// Backward compatibility: If it's a number, it might be an old format hall-recorders.txt
			try {
				double x = std::stod(first_token);
				double y, z, dummy;
				if ((ss >> y >> z) && !(ss >> dummy)) {
					recorders.push_back(std::make_shared<Recorder>(config, (int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_), total_steps__, dir_path));
				}
			} catch (...) {
				// Not a number, not an 'R', ignore line
			}
		}
	}
	file.close();
    return recorders;
}
