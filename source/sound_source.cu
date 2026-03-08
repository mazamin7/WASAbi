#include "sound_source.h"
#include "gaussian_source.h"
#include "simulation.h"
#include "partition.h"
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>


SoundSource::SoundSource(int x, int y, int z, std::string dir_path) :x_(x), y_(y), z_(z)
{
	static int id_generator = 0;
	id_ = id_generator++;
	std::string filename;
	
	filename = dir_path + "/source_data_" + std::to_string(id_) + ".bin";
	source_.open(filename, std::ios::out | std::ios::binary);

}


SoundSource::~SoundSource()
{
}

#include "json.hpp"
using json = nlohmann::json;

std::vector<std::shared_ptr<SoundSource>> SoundSource::ImportSources(std::string path, std::string dir_path)
{
	std::vector<std::shared_ptr<SoundSource>> sources;

    if (std::filesystem::path(path).extension() == ".json") {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "WARNING: Could not open JSON asset file: " << path << std::endl;
            return sources;
        }
        try {
            json j;
            file >> j;
            if (j.contains("sources")) {
                for (auto& s : j["sources"]) {
                    std::string type = "gaussian";
                    if (s.contains("type")) type = s["type"];

                    if (type == "gaussian") {
                        sources.push_back(std::make_shared<GaussianSource>(
                            (int)((double)s["x"] / Simulation::dh_), 
                            (int)((double)s["y"] / Simulation::dh_), 
                            (int)((double)s["z"] / Simulation::dh_),
                            dir_path));
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "ERROR parsing JSON asset file: " << e.what() << std::endl;
        }
        return sources;
    }

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
				sources.push_back(std::make_shared<GaussianSource>((int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_), dir_path));
			}
		}
		else {
			// Backward compatibility: If it's a number, it might be an old format hall-sources.txt
			try {
				double x = std::stod(first_token);
				double y, z, dummy;
				if ((ss >> y >> z) && !(ss >> dummy)) {
					sources.push_back(std::make_shared<GaussianSource>((int)(x / Simulation::dh_), (int)(y / Simulation::dh_), (int)(z / Simulation::dh_), dir_path));
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
	int total_steps = (int)(Simulation::duration_ / Simulation::dt_);
	for (int t = 0; t < total_steps; t++)
	{
		double val = this->SampleValue(t);
		source_.write(reinterpret_cast<const char*>(&val), sizeof(double));
	}
	source_.close();
}
