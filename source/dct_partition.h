#pragma once
#include "partition.h"
#include "dct_volume.h"

class DctPartition : public Partition
{
public:
	DctPartition(int xs, int ys, int zs, int w, int h, int d);
	~DctPartition();

	void Update() override;

	// Spectral components
	DctVolume* pressure_vol_;
	DctVolume* velocity_vol_;
	DctVolume* force_vol_;

private:
	// Precomputed frequency/absorption arrays residing on GPU
	double* d_w0_;
	double* d_alpha_;

	// Pre-calculated physics update coefficients
	// A, B, C: velocity coefficients
	// D, E, F: pressure coefficients
	double* d_coef_A_;
	double* d_coef_B_;
	double* d_coef_C_;
	double* d_coef_D_;
	double* d_coef_E_;
	double* d_coef_F_;

	// Helper function to initialize w0 and alpha arrays on device
	void InitializeConstants();
};
