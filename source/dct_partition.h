#pragma once
#include "partition.h"
#include "dct_volume.h"

class DctPartition : public Partition
{
	double lx2_, ly2_, lz2_;	// actual length ^2
	double* w0_{ nullptr };	// w0
	double* alpha_{ nullptr };	// alpha

	DctVolume pressure_;
	DctVolume velocity_;
	DctVolume residue_;
	DctVolume force_;

	double *next_pressure_modes_{ nullptr };
	double *next_velocity_modes_{ nullptr };

public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d);
	~DctPartition();

	virtual void Update();

	virtual double* get_pressure_field();

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

	virtual std::vector<double> get_xy_forcing_plane(int z);

	std::vector<double> get_xy_force_plane(int z);
	friend class Boundary;
};

