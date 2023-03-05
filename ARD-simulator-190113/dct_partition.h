#pragma once
#include "partition.h"
#include "dct_volume.h"

class DctPartition : public Partition
{
	double lx2_, ly2_, lz2_;	// actual length ^2
	double *cwt_{ nullptr };	// cos(wt)
	double* swt_{ nullptr };	// sin(wt)
	double* alpha_{ nullptr };	// alpha
	double* alpha2_{ nullptr };	// alpha^2
	double* eatm_{ nullptr };	// exp(-alpha*t)
	double* w_omega_{ nullptr };	// w
	double* inv_w_{ nullptr };	// w^-1
	double* inv_w2_{ nullptr };	// w^-2

	bool is_damped_;

	DctVolume pressure_;
	DctVolume force_;
	DctVolume residue_;

	double *prev_modes_{ nullptr };
	double *next_modes_{ nullptr };
	double *velocity_{ nullptr };

public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d, double h_abs);
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

	virtual double check_reset_residues();

	virtual std::vector<double> get_xy_forcing_plane(int z);

	std::vector<double> get_xy_force_plane(int z);
	friend class Boundary;
};

