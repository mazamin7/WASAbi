#pragma once

#include <cufft.h>

class CuDctVolume
{
public:
	CuDctVolume(int w, int h, int d);
	~CuDctVolume();

	void ExecuteDct();
	void ExecuteIdct();
	void reset();
	
	// Data pointers allocated on GPU VRAM
	double* d_values_;
	double* d_modes_;

	int width_;
	int height_;
	int depth_;

private:
	cufftHandle r2c_plan_;
	cufftHandle c2r_plan_;
	
	// Intermediate buffers for cuFFT
	cufftDoubleComplex* d_complex_modes_;
	double* d_extended_;  // 2N-extended real signal for DCT/IDCT
};
