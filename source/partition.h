#pragma once
#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <memory>
//#include "sound_source.h"

class Boundary;
class SoundSource;

class Partition
{
public:
    double dh_;
    double dt_;
    double c0_;
    double air_absorption_alpha1_;
    double air_absorption_alpha2_;

    int x_start_, x_end_;
    int y_start_, y_end_;
    int z_start_, z_end_;
    int width_, height_, depth_;

	struct Info
	{
		int id;
		std::string type;
		int num_sources{ 0 };
		int num_boundaries{ 0 };
	} info_;

	std::vector<std::shared_ptr<SoundSource>> sources_;
	std::vector<std::vector<int>> right_free_borders_;
	std::vector<std::vector<int>> left_free_borders_;
	std::vector<std::vector<int>> top_free_borders_;
	std::vector<std::vector<int>> bottom_free_borders_;
	std::vector<std::vector<int>> front_free_borders_;
	std::vector<std::vector<int>> back_free_borders_;

	bool include_self_terms_{ true };
	bool should_render_{ true };
	bool is_x_pml_{ false };
	bool is_y_pml_{ false };
	bool is_z_pml_{ false };
public:
	static double boundary_absorption_;

	Partition(int xs, int ys, int zs, int w, int h, int d);
	virtual ~Partition();

	virtual void Update() = 0;

    // GPU-resident data pointers (device memory)
    double* d_pressure_;
    double* d_velocity_;
    double* d_force_;
    double* d_residue_;
    double* d_max_p_; // for dynamic scaling reduction

    // Array of source metadata
    int num_sources_ = 0;
    int* d_source_indices_ = nullptr; 
    double* d_source_values_ = nullptr;

    cudaStream_t stream_;

	virtual double* get_pressure_field();
	virtual std::vector<double> get_xy_plane(int z);
	virtual std::vector<double> get_yz_plane(int x);
	virtual std::vector<double> get_xz_plane(int y);

    virtual double get_pressure(int x, int y, int z);
    virtual void set_pressure(int x, int y, int z, double v);
    virtual void add_to_pressure(int x, int y, int z, double v);

    virtual double get_velocity(int x, int y, int z);
    virtual void set_velocity(int x, int y, int z, double v);
    virtual void add_to_velocity(int x, int y, int z, double v);

    virtual double get_residue(int x, int y, int z);
    virtual void set_residue(int x, int y, int z, double v);
    virtual void add_to_residue(int x, int y, int z, double v);

    virtual double get_force(int x, int y, int z);
    virtual void set_force(int x, int y, int z, double v);

	virtual void reset_forces();
	virtual void reset_residues();

	virtual void PostMerge();
	virtual std::vector<double> get_xy_forcing_plane(int z);

	void AddBoundary(std::shared_ptr<Boundary> boundary);
	void AddBoundary(Boundary* boundary);
	void AddSource(std::shared_ptr<SoundSource> source);
	static std::vector<std::shared_ptr<Partition>> ImportPartitions(std::string path);
	void Info();

	virtual void ComputeSourceForcingTerms(double t);
	virtual double GetMaxAbsolutePressure();

    void RenderToBuffer(uint32_t* d_pixels, int plane_type, int coord, int screen_width, int screen_height, int x_offset, int y_offset, float v_coef);

    std::vector<Boundary*> cu_boundaries_;

	friend class Boundary;
	friend class Simulation;
	friend class Tools;
	friend class PmlPartition;
	friend class Recorder;
};

