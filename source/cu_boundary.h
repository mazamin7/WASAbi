#pragma once
#include <memory>

#include "boundary.h"
#include "cu_partition.h"

class CuPartition;

class CuBoundary : public Boundary
{
public:
	CuBoundary(BoundaryType type, double absorp, std::shared_ptr<CuPartition> a, std::shared_ptr<CuPartition> b,
		int xs, int xe, int ys, int ye, int zs, int ze);
	virtual ~CuBoundary();

	virtual void ComputeResidues() override;

	double* d_coefs_; // Device pointer for coefficients matrix

	static std::shared_ptr<CuBoundary> FindBoundary(std::shared_ptr<CuPartition> a, std::shared_ptr<CuPartition> b, double absorp = 1.0);
};
