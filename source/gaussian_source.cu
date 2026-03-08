#define _USE_MATH_DEFINES
#include <cmath>
#include "gaussian_source.h"
#include "simulation.h"


GaussianSource::GaussianSource(int x, int y, int z, std::string dir_path) :SoundSource(x, y, z, dir_path)
{
}


GaussianSource::~GaussianSource()
{
}

double GaussianSource::SampleValue(double t)
{
	double arg = pow(M_PI * ((2 * (Simulation::c0_*Simulation::dt_ / Simulation::dh_) * t) / 6 - 2.0), 2);
	return -1e9 * exp(-arg);
}
