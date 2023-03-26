#pragma once
#include "partition.h"
#include "dct_volume.h"

class DctPartition : public Partition
{
	double lx2_, ly2_, lz2_;	// actual length ^2
	double *cwt_{ nullptr };	// cos(wt)
	double* swt_{ nullptr };	// sin(wt)
	double* w_omega_{ nullptr };	// w
	double* w2_{ nullptr };		// w^2
	double* inv_w_{ nullptr };	// w^-1
	double* inv_w2_{ nullptr };	// w^-2
	double alpha_;	// alpha
	double alpha2_;	// alpha^2
	double eatm_;	// exp(-alpha*t)

	bool is_damped_;

	DctVolume pressure_;
	DctVolume force_;
	DctVolume residue_;

	double *prev_modes_{ nullptr };
	double *next_modes_{ nullptr };
	double *velocity_{ nullptr };

public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d, double alpha_abs);
	~DctPartition();

	virtual void Update();

	virtual double* get_pressure_field();
	virtual double get_pressure(int x, int y, int z);
	virtual void set_pressure(int x, int y, int z, double v);
	virtual double get_residue(int x, int y, int z);
	virtual void set_residue(int x, int y, int z, double v);
	virtual double get_force(int x, int y, int z);
	virtual void set_force(int x, int y, int z, double f);

	virtual void reset_forces();
	virtual void reset_residues();

	virtual void add_to_pressure(int x, int y, int z, double v);

	virtual std::vector<double> get_xy_forcing_plane(int z);

	std::vector<double> get_xy_force_plane(int z);
	friend class Boundary;
};

