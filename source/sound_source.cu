#include "sound_source.h"
#include "gaussian_source.h"
#include "simulation.h"
#include "partition.h"
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>


SoundSource::SoundSource(int x, int y, int z) :x_(x), y_(y), z_(z)
{
	static int id_generator = 0;
	id_ = id_generator++;
	std::string filename;
	std::string dir_name = std::to_string(Simulation::dh_) + "_" + std::to_string(Partition::boundary_absorption_) + "_" + std::to_string(Simulation::air_absorption_alpha1_) + "_" + std::to_string(Simulation::air_absorption_alpha2_);
	
	std::string dir_path = "./output/" + dir_name;
	std::filesystem::create_directories(dir_path);
	
	filename = dir_path + "/source_" + std::to_string(id_) + ".txt";
	source_.open(filename, std::ios::out);

}


SoundSource::~SoundSource()
{
}

std::vector<std::shared_ptr<SoundSource>> SoundSource::ImportSources(std::string path)
{
	std::vector<std::shared_ptr<SoundSource>> sources;
	std::ifstream file(path);
	if (!file.is_open())
	{
		std::cerr << "WARNING: Could not open unified asset file: " << path << std::endl;
		return sources;
	}

	std::string line;
	while (std::getline(file, line))
	{
		if (line.empty()) continue;

		std::stringstream ss(line);
		std::string first_token;
		ss >> first_token;

		if (first_token == "S" || first_token == "s") {
			// It's explicitly a source
			double x, y, z;
			if (ss >> x >> y >> z) {
				sources.push_back(std::make_shared<GaussianSource>((int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_)));
			}
		}
		else {
			// Backward compatibility: If it's a number, it might be an old format hall-sources.txt
			try {
				double x = std::stod(first_token);
				double y, z, dummy;
				if ((ss >> y >> z) && !(ss >> dummy)) {
					sources.push_back(std::make_shared<GaussianSource>((int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_)));
				}
			} catch (...) {
				// Not a number, not an 'S', ignore line
			}
		}
	}
	file.close();
	for (auto source : sources)
	{
		source->RecordSource();
	}
	return sources;
}

void SoundSource::RecordSource()
{
	for (int t = 0; t < Simulation::duration_ / Simulation::dt_; t++)
	{
		source_ << this->SampleValue(t) << std::endl;
	}
	source_.close();
}
