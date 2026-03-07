#pragma once
#include "cu_partition.h"
#include "cu_dct_volume.h"

class CuDctPartition : public CuPartition
{
public:
	CuDctPartition(int xs, int ys, int zs, int w, int h, int d);
	~CuDctPartition();

	void Update() override;

	// Spectral components
	CuDctVolume* pressure_vol_;
	CuDctVolume* velocity_vol_;
	CuDctVolume* force_vol_;

private:
	// Precomputed frequency/absorption arrays residing on GPU
	double* d_w0_;
	double* d_alpha_;

	// Helper function to initialize w0 and alpha arrays on device
	void InitializeConstants();
};
