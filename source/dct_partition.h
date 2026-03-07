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

	// Helper function to initialize w0 and alpha arrays on device
	void InitializeConstants();
};
