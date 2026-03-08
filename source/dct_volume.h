#pragma once

#include <cufft.h>

class DctVolume
{
public:
	DctVolume(int w, int h, int d, double* shared_ext = nullptr, cufftDoubleComplex* shared_complex = nullptr);
	~DctVolume();

	void ExecuteDct(cudaStream_t stream = 0);
	void ExecuteIdct(cudaStream_t stream = 0);
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
	
	bool owns_extended_;
	bool owns_complex_;
	cufftDoubleComplex* d_complex_modes_;
	double* d_extended_;  // 2N-extended real signal for DCT/IDCT
};
