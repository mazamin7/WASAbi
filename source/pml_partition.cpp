#include "pml_partition.h"
#include "simulation.h"
#include <omp.h>


int PmlPartition::GetIndex(int x, int y, int z)
{
	if (x < 0 || x >= width_)
		return width_ * height_ * depth_;
	if (y < 0 || y >= height_)
		return width_ * height_ * depth_;
	if (z < 0 || z >= depth_)
		return width_ * height_ * depth_;
	return z * height_ * width_ + y * width_ + x;
}

PmlPartition::PmlPartition(std::shared_ptr<Partition> neighbor_part, PmlType type, int xs, int ys, int zs, int w, int h, int d)
	: Partition(xs, ys, zs, w, h, d), type_(type), neighbor_part_(neighbor_part)
{
	include_self_terms_ = true;
	should_render_ = false;
	info_.type = "PML";

	if (type_ == P_LEFT || type_ == P_RIGHT) is_x_pml_ = true;
	if (type_ == P_TOP || type_ == P_BOTTOM) is_y_pml_ = true;
	if (type_ == P_FRONT || type_ == P_BACK) is_z_pml_ = true;

	thickness_ = Simulation::n_pml_layers_ * dh_;
	zeta_ = Simulation::c0_ / thickness_ * log10(1 / R_);

	int size = width_ * height_ * depth_ + 1;

	/*p_ = (double *)malloc(size * sizeof(double));
	p_new_ = (double *)malloc(size * sizeof(double));

	v_ = (double*)malloc(size * sizeof(double));
	v_new_ = (double*)malloc(size * sizeof(double));

	psi_ = (double*)malloc(size * sizeof(double));
	psi_new_ = (double*)malloc(size * sizeof(double));

	phi_x_ = (double *)malloc(size * sizeof(double));
	phi_x_new_ = (double *)malloc(size * sizeof(double));
	phi_y_ = (double *)malloc(size * sizeof(double));
	phi_y_new_ = (double *)malloc(size * sizeof(double));
	phi_z_ = (double *)malloc(size * sizeof(double));
	phi_z_new_ = (double *)malloc(size * sizeof(double));

	residue_ = (double*)malloc(size * sizeof(double));
	force_ = (double *)malloc(size * sizeof(double));

	memset((void *)p_, 0, size * sizeof(double));
	memset((void *)p_new_, 0, size * sizeof(double));
	memset((void*)v_, 0, size * sizeof(double));
	memset((void*)v_new_, 0, size * sizeof(double));
	memset((void*)psi_, 0, size * sizeof(double));
	memset((void*)psi_new_, 0, size * sizeof(double));
	memset((void *)phi_x_, 0, size * sizeof(double));
	memset((void *)phi_x_new_, 0, size * sizeof(double));
	memset((void *)phi_y_, 0, size * sizeof(double));
	memset((void *)phi_y_new_, 0, size * sizeof(double));
	memset((void *)phi_z_, 0, size * sizeof(double));
	memset((void *)phi_z_new_, 0, size * sizeof(double));
	memset((void *)force_, 0, size * sizeof(double));
	memset((void *)residue_, 0, size * sizeof(double));*/

	p_ = std::make_unique<double[]>(size);
	p_new_ = std::make_unique<double[]>(size);
	v_ = std::make_unique<double[]>(size);
	v_new_ = std::make_unique<double[]>(size);
	psi_ = std::make_unique<double[]>(size);
	psi_new_ = std::make_unique<double[]>(size);
	phi_x_ = std::make_unique<double[]>(size);
	phi_x_new_ = std::make_unique<double[]>(size);
	phi_y_ = std::make_unique<double[]>(size);
	phi_y_new_ = std::make_unique<double[]>(size);
	phi_z_ = std::make_unique<double[]>(size);
	phi_z_new_ = std::make_unique<double[]>(size);
	residue_ = std::make_unique<double[]>(size);
	force_ = std::make_unique<double[]>(size);
	zetax_ = std::make_unique<double[]>(size);
	zetay_ = std::make_unique<double[]>(size);
	zetaz_ = std::make_unique<double[]>(size);

	std::memset(p_.get(), 0, size * sizeof(double));
	std::memset(p_new_.get(), 0, size * sizeof(double));
	std::memset(v_.get(), 0, size * sizeof(double));
	std::memset(v_new_.get(), 0, size * sizeof(double));
	std::memset(psi_.get(), 0, size * sizeof(double));
	std::memset(psi_new_.get(), 0, size * sizeof(double));
	std::memset(phi_x_.get(), 0, size * sizeof(double));
	std::memset(phi_x_new_.get(), 0, size * sizeof(double));
	std::memset(phi_y_.get(), 0, size * sizeof(double));
	std::memset(phi_y_new_.get(), 0, size * sizeof(double));
	std::memset(phi_z_.get(), 0, size * sizeof(double));
	std::memset(phi_z_new_.get(), 0, size * sizeof(double));
	std::memset(force_.get(), 0, size * sizeof(double));
	std::memset(residue_.get(), 0, size * sizeof(double));

#pragma omp parallel for collapse(2)
	for (int k = 0; k < depth_; k++)
	{
		for (int j = 0; j < height_; j++)
		{
			for (int i = 0; i < width_; i++)
			{
				switch (type)
				{
				case PmlPartition::P_BACK:
					zetaz_[GetIndex(i, j, k)] = zeta_ * ((k + 1) * dh_ / thickness_ - sin(2 * M_PI * (k + 1) * dh_ / thickness_) / 2 / M_PI);
					break;
				case PmlPartition::P_FRONT:
					zetaz_[GetIndex(i, j, k)] = zeta_ * ((depth_ - k) * dh_ / thickness_ - sin(2 * M_PI * (depth_ - k) * dh_ / thickness_) / 2 / M_PI);
					break;
				case PmlPartition::P_BOTTOM:
					zetay_[GetIndex(i, j, k)] = zeta_ * ((j + 1) * dh_ / thickness_ - sin(2 * M_PI * (j + 1) * dh_ / thickness_) / 2 / M_PI);
					break;
				case PmlPartition::P_TOP:
					zetay_[GetIndex(i, j, k)] = zeta_ * ((height_ - j) * dh_ / thickness_ - sin(2 * M_PI * (height_ - j) * dh_ / thickness_) / 2 / M_PI);
					break;
				case PmlPartition::P_RIGHT:
					// zeta_x = zeta * ( (|x - a|)/L - (sin(2 Pi (|x - a|)/(L))) / (2 Pi) )
					zetax_[GetIndex(i, j, k)] = zeta_ * ((i + 1) * dh_ / thickness_ - sin(2 * M_PI * (i + 1) * dh_ / thickness_) / 2 / M_PI);
					break;
				case PmlPartition::P_LEFT:
					zetax_[GetIndex(i, j, k)] = zeta_ * ((width_ - i) * dh_ / thickness_ - sin(2 * M_PI * (width_ - i) * dh_ / thickness_) / 2 / M_PI);
					break;
				default:
					break;
				}
			}
		}
	}
}


PmlPartition::~PmlPartition()
{
	/*free(p_);
	free(v_);
	free(v_new_);
	free(psi_);
	free(psi_new_);
	free(phi_x_);
	free(phi_x_new_);
	free(phi_y_);
	free(phi_y_new_);
	free(phi_z_);
	free(phi_z_new_);
	free(psi_);
	free(psi_new_);
	free(force_);
	free(residue_);
	free(zetax_);
	free(zetay_);
	free(zetaz_);
	free(p_new_);*/
}

void PmlPartition::Update()
{
	int width = width_;
	int height = height_;
	int depth = depth_;
	auto type = type_;
	auto thickness = thickness_;
	auto dh = Simulation::dh_;
	auto dt = Simulation::dt_;
	auto c0 = Simulation::c0_;
	auto zeta = zeta_;

	double coefs[7];
	double amp = 0.0;

	double temp_coefs[] = { 2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0 };
	std::copy(std::begin(temp_coefs), std::end(temp_coefs), std::begin(coefs));
	amp = 180.0;

	double fourthCoefs[] = { 1.0, -8.0, 0.0, 8.0, -1.0 };

#pragma omp parallel for collapse(2)
	for (int k = 0; k < depth; k++)
	{
		//#pragma omp parallel for
		for (int j = 0; j < height; j++)
		{
			//#pragma omp parallel for
			for (int i = 0; i < width; i++)
			{
				p_new_[GetIndex(i, j, k)] = p_[GetIndex(i, j, k)] + dt_ * v_[GetIndex(i, j, k)];
			}
		}
	}

#pragma omp parallel for collapse(2)
	for (int k = 0; k < depth; k++)
	{
		//#pragma omp parallel for
		for (int j = 0; j < height; j++)
		{
			//#pragma omp parallel for
			for (int i = 0; i < width; i++)
			{
				psi_new_[GetIndex(i, j, k)] = psi_[GetIndex(i, j, k)] + dt * p_new_[GetIndex(i, j, k)];

				// --------------------------------------------------------

				double d2udx2 = 0.0;
				double d2udy2 = 0.0;
				double d2udz2 = 0.0;

				for (int m = 0; m < 7; m++)
				{
					d2udx2 += coefs[m] * p_new_[GetIndex(i + m - 3, j, k)];
					d2udy2 += coefs[m] * p_new_[GetIndex(i, j + m - 3, k)];
					d2udz2 += coefs[m] * p_new_[GetIndex(i, j, k + m - 3)];
				}

				d2udx2 /= (amp * dh * dh);
				d2udy2 /= (amp * dh * dh);
				d2udz2 /= (amp * dh * dh);


				// v_{t}
				double term1 = v_[GetIndex(i, j, k)];

				// c_0^2 Delta p
				double term3 = c0 * c0 * (d2udx2 + d2udy2 + d2udz2);	// c^2*(d2udx2+d2udy2+d2udz2)

				// -(Zx + Zy + Zz) p_t
				double term4 = -(zetax_[GetIndex(i, j, k)] + zetay_[GetIndex(i, j, k)] + zetaz_[GetIndex(i, j, k)]) * v_[GetIndex(i, j, k)] / 2;

				// -(Zy Zz + Zx Zz + Zy Zz) p
				double term5 = -(zetax_[GetIndex(i, j, k)] * zetay_[GetIndex(i, j, k)] + zetay_[GetIndex(i, j, k)] * zetaz_[GetIndex(i, j, k)] + zetax_[GetIndex(i, j, k)] * zetaz_[GetIndex(i, j, k)]) * p_new_[GetIndex(i, j, k)];


				double dphidx = 0.0;
				double dphidy = 0.0;
				double dphidz = 0.0;

				for (int m = 0; m < 5; m++)
				{
					dphidx += fourthCoefs[m] * phi_x_[GetIndex(i + m - 2, j, k)];
					dphidy += fourthCoefs[m] * phi_y_[GetIndex(i, j + m - 2, k)];
					dphidz += fourthCoefs[m] * phi_z_[GetIndex(i, j, k + m - 2)];
				}

				dphidx /= (12.0 * dh);
				dphidy /= (12.0 * dh);
				dphidz /= (12.0 * dh);


				// Nabla phi
				double term6 = dphidx + dphidy + dphidz;

				// -Zx Zy Zz psi
				double term7 = -zetax_[GetIndex(i, j, k)] * zetay_[GetIndex(i, j, k)] * zetaz_[GetIndex(i, j, k)] * (psi_new_[GetIndex(i, j, k)] + psi_[GetIndex(i, j, k)]) * 0.5;


				// v_{t} = c_0^2 Delta p - (Zx + Zy + Zz) p_t - (Zy Zz + Zx Zz + Zy Zz) p + Nabla phi - Zx Zy Zz psi
				v_new_[GetIndex(i, j, k)] = term1 + dt * (term3 + term4 + term5 + term6 + term7 + force_[GetIndex(i, j, k)]);
				v_new_[GetIndex(i, j, k)] = v_new_[GetIndex(i, j, k)] / (1 + dt / 2 * (zetax_[GetIndex(i, j, k)] + zetay_[GetIndex(i, j, k)] + zetaz_[GetIndex(i, j, k)]));

				// --------------------------------------------------------

				double dudx = 0.0;
				double dudy = 0.0;
				double dudz = 0.0;

				for (int m = 0; m < 5; m++)
				{
					dudx += fourthCoefs[m] * p_new_[GetIndex(i + m - 2, j, k)];
					dudy += fourthCoefs[m] * p_new_[GetIndex(i, j + m - 2, k)];
					dudz += fourthCoefs[m] * p_new_[GetIndex(i, j, k + m - 2)];
				}

				dudx /= (12.0 * dh);
				dudy /= (12.0 * dh);
				dudz /= (12.0 * dh);


				double dpsidx = 0.0;
				double dpsidy = 0.0;
				double dpsidz = 0.0;

				for (int m = 0; m < 5; m++)
				{
					dpsidx += fourthCoefs[m] * p_new_[GetIndex(i + m - 2, j, k)];
					dpsidy += fourthCoefs[m] * p_new_[GetIndex(i, j + m - 2, k)];
					dpsidz += fourthCoefs[m] * p_new_[GetIndex(i, j, k + m - 2)];
				}

				dpsidx /= (12.0 * dh);
				dpsidy /= (12.0 * dh);
				dpsidz /= (12.0 * dh);


				// phi_t
				double phi_x_term_1 = phi_x_[GetIndex(i, j, k)];
				double phi_y_term_1 = phi_y_[GetIndex(i, j, k)];
				double phi_z_term_1 = phi_z_[GetIndex(i, j, k)];

				// Gamma_11 phi
				double phi_x_term_2 = -zetax_[GetIndex(i, j, k)] * phi_x_[GetIndex(i, j, k)];
				double phi_y_term_2 = -zetay_[GetIndex(i, j, k)] * phi_y_[GetIndex(i, j, k)];
				double phi_z_term_2 = -zetaz_[GetIndex(i, j, k)] * phi_z_[GetIndex(i, j, k)];

				// c_0^2 Gamma_2 Nabla p
				double phi_x_term_3 = c0 * c0 * (zetay_[GetIndex(i, j, k)] + zetaz_[GetIndex(i, j, k)] - zetax_[GetIndex(i, j, k)]) * dudx;
				double phi_y_term_3 = c0 * c0 * (zetax_[GetIndex(i, j, k)] + zetaz_[GetIndex(i, j, k)] - zetay_[GetIndex(i, j, k)]) * dudy;
				double phi_z_term_3 = c0 * c0 * (zetax_[GetIndex(i, j, k)] + zetay_[GetIndex(i, j, k)] - zetaz_[GetIndex(i, j, k)]) * dudz;

				// c_0^2 Gamma_3 Nabla psi
				double phi_x_term_4 = c0 * c0 * zetay_[GetIndex(i, j, k)] * zetaz_[GetIndex(i, j, k)] * dpsidx;
				double phi_y_term_4 = c0 * c0 * zetax_[GetIndex(i, j, k)] * zetaz_[GetIndex(i, j, k)] * dpsidy;
				double phi_z_term_4 = c0 * c0 * zetax_[GetIndex(i, j, k)] * zetay_[GetIndex(i, j, k)] * dpsidz;


				// phi_t = Gamma_1 phi + c_0^2 Gamma_2 Nabla p + c_0^2 Gamma_3 Nabla psi
				phi_x_new_[GetIndex(i, j, k)] = phi_x_term_1 + dt * (phi_x_term_2 + phi_x_term_3 + phi_x_term_4);
				phi_y_new_[GetIndex(i, j, k)] = phi_y_term_1 + dt * (phi_y_term_2 + phi_y_term_3 + phi_y_term_4);
				phi_z_new_[GetIndex(i, j, k)] = phi_z_term_1 + dt * (phi_z_term_2 + phi_z_term_3 + phi_z_term_4);
			}
		}
	}

	std::swap(v_, v_new_);
	std::swap(p_, p_new_);
	std::swap(psi_, psi_new_);
	std::swap(phi_x_, phi_x_new_);
	std::swap(phi_y_, phi_y_new_);
	std::swap(phi_z_, phi_z_new_);
}

double* PmlPartition::get_pressure_field()
{
	return p_.get();
}

double PmlPartition::get_pressure(int x, int y, int z)
{
	return p_[GetIndex(x, y, z)];
}

void PmlPartition::set_pressure(int x, int y, int z, double v)
{
	p_[GetIndex(x, y, z)] = v;
}

void PmlPartition::add_to_pressure(int x, int y, int z, double v)
{
	p_[GetIndex(x, y, z)] = p_[GetIndex(x, y, z)] + v;
}

double PmlPartition::get_velocity(int x, int y, int z)
{
	return v_[GetIndex(x, y, z)];
}

void PmlPartition::set_velocity(int x, int y, int z, double v)
{
	v_[GetIndex(x, y, z)] = v;
}

void PmlPartition::add_to_velocity(int x, int y, int z, double v)
{
	v_[GetIndex(x, y, z)] = v_[GetIndex(x, y, z)] + v;
}

double PmlPartition::get_residue(int x, int y, int z)
{
	return residue_[GetIndex(x, y, z)];
}

void PmlPartition::set_residue(int x, int y, int z, double v)
{
	residue_[GetIndex(x, y, z)] = v;
}

void PmlPartition::add_to_residue(int x, int y, int z, double v)
{
	// #pragma omp atomic
	residue_[GetIndex(x, y, z)] = residue_[GetIndex(x, y, z)] + v;
}

double PmlPartition::get_force(int x, int y, int z)
{
	return force_[GetIndex(x, y, z)];
}

void PmlPartition::set_force(int x, int y, int z, double f)
{
	force_[GetIndex(x, y, z)] = f;
}

void PmlPartition::reset_forces()
{
	int width = width_;
	int height = height_;
	int depth = depth_;

	memset((void*)force_.get(), 0, width * height * depth * sizeof(double));
}

void PmlPartition::reset_residues()
{
	int width = width_;
	int height = height_;
	int depth = depth_;

	memset((void*)residue_.get(), 0, width * height * depth * sizeof(double));
}
