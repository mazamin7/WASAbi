#pragma once
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include "partition.h"
#include "config_loader.h"

class Recorder
{
	int id_;
	int x_, y_, z_;
	int total_steps_;
	Config config_;

public:
	int x() const { return x_; }
	int y() const { return y_; }
	int z() const { return z_; }

private:
	std::shared_ptr<Partition> part_;
	std::vector<std::shared_ptr<Partition>> partitions_;
	
	std::fstream output_;
	std::fstream response_;
	std::string output_path_;
	std::string response_path_;
	std::string dir_path_;

	double* field_buffer_ = nullptr;
	int gs_x_ = 0, gs_y_ = 0, gs_z_ = 0;
	int go_x_ = 0, go_y_ = 0, go_z_ = 0;

public:
	Recorder(const Config& config, int x, int y, int z, int total_steps__, std::string dir_path);
	~Recorder();

	void FindPartition(std::vector<std::shared_ptr<Partition>> partitions);
	void RecordField(int time_step = 0);
	void RecordResponse(int time_step = 0);

	static std::vector<std::shared_ptr<Recorder>> ImportRecorders(const Config& config, std::string path, int total_steps__, std::string dir_path);

};

