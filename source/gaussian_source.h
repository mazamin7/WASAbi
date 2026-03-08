#pragma once
#include "sound_source.h"

class GaussianSource :public SoundSource
{
public:
	GaussianSource(int x, int y,int z, std::string dir_path);
	~GaussianSource();

	virtual double SampleValue(double t);
};

