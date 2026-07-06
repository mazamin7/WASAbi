#pragma once
#include "partition.h"
#include "dct_volume.h"

class DctPartition : public Partition
{
public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d);
	~DctPartition();

	void Update() override;

	void HalfKick();
	void Drift();
	void MergeResiduesIntoForce();

	// Spectral components
	DctVolume* pressure_vol_;
	DctVolume* velocity_vol_;
	DctVolume* force_vol_;

private:
	// Precomputed frequency/absorption arrays residing on GPU
	double* d_w0_;
	double* d_alpha_;
	double* d_shared_ext_;
	cufftDoubleComplex* d_shared_complex_;

	// Pre-calculated physics update coefficients
	// A, B, C: velocity coefficients
	// D, E, F: pressure coefficients
	// Pre-calculated exact drift matrix S_m
	double* d_S11_;
	double* d_S12_;
	double* d_S21_;
	double* d_S22_;

	// Helper function to initialize w0 and alpha arrays on device
	void InitializeConstants();
};
