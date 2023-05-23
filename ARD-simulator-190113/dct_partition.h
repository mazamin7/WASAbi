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

	bool second_order_; // true -> second order, false -> first order

	DctVolume pressure_;
	DctVolume velocity_;
	DctVolume residue_;
	DctVolume force_;

	double *prev_modes_{ nullptr };
	double *next_modes_{ nullptr };
	double *velocity_modal_{ nullptr };

public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d, double alpha_abs);
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
	virtual void add_to_force(int x, int y, int z, double v);

	virtual void reset_forces();
	virtual void reset_residues();

	virtual std::vector<double> get_xy_forcing_plane(int z);

	std::vector<double> get_xy_force_plane(int z);
	friend class Boundary;
};

