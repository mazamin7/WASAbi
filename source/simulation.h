#pragma once
#include <vector>
#include <memory>
#include <string>
#include <SDL.h>

class Partition;
class Boundary;
class SoundSource;

class Simulation
{
	std::vector<std::shared_ptr<Partition>> partitions_;
	std::vector<std::shared_ptr<Partition>> dct_partitions_;
	std::vector<std::shared_ptr<Partition>> pml_partitions_;
	std::vector<std::shared_ptr<Boundary>> boundaries_;
	std::vector<std::shared_ptr<SoundSource>> sources_;

	int x_start_, x_end_;
	int y_start_, y_end_;
	int z_start_, z_end_;

	int size_x_, size_y_, size_z_;

	bool ready_;
	std::vector<Uint32> pixels_;
	uint32_t* d_pixels_ = nullptr;
	SDL_PixelFormat* sdl_fmt_;
	double max_p_{ 0.0 };

	struct Info
	{
		size_t num_partitions{ 0 };
		size_t num_dct_partitions{ 0 };
		size_t num_pml_partitions{ 0 };
		size_t num_sources{ 0 };
		size_t num_boundaries{ 0 };
		std::vector<std::vector<char>> model_map;
	} info_;

public:

	static double duration_;

	static double dh_;
	static double dt_;
	static double c0_;
	static double air_absorption_alpha1_;
	static double air_absorption_alpha2_;
	static int n_pml_layers_;
	static int viz_skip_;
	static float max_viz_gain_;

	int time_step_{ 0 };
	int panel_w_{ 0 }, panel_h_{ 0 }; // size of each individual panel
	float v_coef_{ 0.1f };           // dynamic scaling factor

	Simulation(std::vector<std::shared_ptr<Partition>> &partitions, std::vector<std::shared_ptr<SoundSource>> &sources);
	~Simulation();

	int Update();
	double GetLastMaxPressure() { return max_p_; }

	void Info();
	//FindBoundaries();

	int size_x() { return size_x_; }
	int size_y() { return size_y_; }
	int size_z() { return size_z_; }
	// Full 3-panel render buffer dimensions (3 panels wide)
	int render_w() { return panel_w_ * 3; }
	int render_h() { return panel_h_; }
	bool ready() { return ready_; }
	decltype(pixels_)& pixels() { return pixels_; }
};

