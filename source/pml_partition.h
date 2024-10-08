#pragma once
#include "partition.h"
#include <memory>

class PmlPartition : public Partition
{
	std::shared_ptr<Partition> neighbor_part_;

	double R_{ 1.0E-5 };
	double zeta_;
	double thickness_;

	std::unique_ptr<double[]> p_;
	std::unique_ptr<double[]> p_new_;

	std::unique_ptr<double[]> v_;
	std::unique_ptr<double[]> v_new_;

	std::unique_ptr<double[]> psi_;
	std::unique_ptr<double[]> psi_new_;

	std::unique_ptr<double[]> phi_x_;
	std::unique_ptr<double[]> phi_x_new_;
	std::unique_ptr<double[]> phi_y_;
	std::unique_ptr<double[]> phi_y_new_;
	std::unique_ptr<double[]> phi_z_;
	std::unique_ptr<double[]> phi_z_new_;

	std::unique_ptr<double[]> zetax_;
	std::unique_ptr<double[]> zetay_;
	std::unique_ptr<double[]> zetaz_;

	std::unique_ptr<double[]> residue_;
	std::unique_ptr<double[]> force_;

	int GetIndex(int x, int y, int z);

public:
	enum PmlType {
		P_LEFT,
		P_RIGHT,
		P_TOP,
		P_BOTTOM,
		P_FRONT,
		P_BACK
	}type_;

	PmlPartition(std::shared_ptr<Partition> neighbor_part, PmlType type, int xs, int ys, int zs, int w, int h, int d);
	~PmlPartition();

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

	virtual void reset_forces();
	virtual void reset_residues();
};

//#pragma once
//#include "partition.h"
//
//class PmlPartition : public Partition
//{
//	std::shared_ptr<Partition> neighbor_part_;
//
//	double R_{1.0E-5};
//	double zeta_;
//	double thickness_;
//
//	double* p_{ nullptr };
//	double* p_new_{ nullptr };
//
//	double* v_{ nullptr };
//	double* v_new_{ nullptr };
//
//	double* psi_;
//	double* psi_new_;
//
//	double* phi_x_;
//	double* phi_x_new_;
//	double* phi_y_;
//	double* phi_y_new_;
//	double* phi_z_;
//	double* phi_z_new_;
//
//	double* zetax_;
//	double* zetay_;
//	double* zetaz_;
//
//	double* residue_{ nullptr };
//	double* force_;
//
//	int GetIndex(int x, int y, int z);
//
//public:
//	enum PmlType {
//		P_LEFT,
//		P_RIGHT,
//		P_TOP,
//		P_BOTTOM,
//		P_FRONT,
//		P_BACK
//	}type_;
//
//	PmlPartition(std::shared_ptr<Partition> neighbor_part, PmlType type, int xs, int ys, int zs, int w, int h, int d);
//	~PmlPartition();
//
//	virtual void Update();
//
//	virtual double* get_pressure_field();
//
//	virtual double get_pressure(int x, int y, int z);
//	virtual void set_pressure(int x, int y, int z, double v);
//	virtual void add_to_pressure(int x, int y, int z, double v);
//
//	virtual double get_velocity(int x, int y, int z);
//	virtual void set_velocity(int x, int y, int z, double v);
//	virtual void add_to_velocity(int x, int y, int z, double v);
//
//	virtual double get_residue(int x, int y, int z);
//	virtual void set_residue(int x, int y, int z, double v);
//	virtual void add_to_residue(int x, int y, int z, double v);
//
//	virtual double get_force(int x, int y, int z);
//	virtual void set_force(int x, int y, int z, double v);
//
//	virtual void reset_forces();
//	virtual void reset_residues();
//};
//



