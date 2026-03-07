#include "simulation.h"
#include "partition.h"
#include "pml_partition.h"
#include <cuda_runtime.h>
#include "boundary.h"

#include "tools.h"
#include "sound_source.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <assert.h>
#include <sstream>


Simulation::Simulation(std::vector<std::shared_ptr<Partition>> &partitions, std::vector<std::shared_ptr<SoundSource>> &sources)
	: partitions_(partitions), sources_(sources)
{
	// Find all shared boundaries of partitions
	for (int i = 0; i < partitions_.size(); i++)
	{
		auto part_a = partitions_[i];
		for (int j = i + 1; j < partitions_.size(); j++)
		{
			auto part_b = partitions_[j];
			auto boundary = Boundary::FindBoundary(part_a, part_b);
			if (boundary)
			{
				boundaries_.push_back(boundary);
				part_a->AddBoundary(boundary);
				part_b->AddBoundary(boundary);
			}
		}
	}

	// Add sources to corresponding partition
	for (auto partition : partitions_)
	{
		for (auto source : sources_)
		{
			if (source->x_ >= partition->x_start_ && source->x_ < partition->x_end_ &&
				source->y_ >= partition->y_start_ && source->y_ < partition->y_end_ &&
				source->z_ >= partition->z_start_ && source->z_ < partition->z_end_)
			{
				partition->AddSource(source);
			}
		}
	}

	info_.num_dct_partitions = partitions_.size();
	info_.num_boundaries = boundaries_.size();
	info_.num_sources = sources_.size();

	dct_partitions_ = partitions_;


	// Find and create PML partitions.
	for (int cnt = 0; cnt < info_.num_dct_partitions; cnt++)
	{
		auto partition = partitions_[cnt];
		int start;
		int end;
		bool started;

		// Find left PML.
		start = 0; end = 0; started = false;
		for (int i = 0; i < partition->height_; i++)
		{
			if (started)
			{
				end++;
				if (partition->left_free_borders_[0][i] && i != partition->height_ - 1)
				{
					continue;
				}
				else if (start != end)
				{
					if (i != partition->height_ - 1) end--;
					auto pml = std::make_shared<PmlPartition>(
						partition,
						PmlPartition::P_LEFT,
						partition->x_start_ - Simulation::n_pml_layers_,
						partition->y_start_ + start,
						partition->z_start_,
						Simulation::n_pml_layers_,
						end - start + 1,
						partition->depth_);
					partitions_.push_back(pml);
					auto boundary = Boundary::FindBoundary(pml, partition, partition->boundary_absorption_);
					boundaries_.push_back(boundary);
					info_.num_pml_partitions++;
					started = false;
				}
			}
			else if (partition->left_free_borders_[0][i])
			{
				start = i;
				end = i;
				started = true;
			}
		}

		// Find right PML.
		start = 0; end = 0; started = false;
		for (int i = 0; i < partition->height_; i++)
		{
			if (started)
			{
				end++;
				if (partition->right_free_borders_[0][i] && i != partition->height_ - 1)
				{
					continue;
				}
				else if (start != end)
				{
					if (i != partition->height_ - 1) end--;
					auto pml = std::make_shared<PmlPartition>(
						partition,
						PmlPartition::P_RIGHT,
						partition->x_end_,
						partition->y_start_ + start,
						partition->z_start_,
						Simulation::n_pml_layers_,
						end - start + 1,
						partition->depth_);
					partitions_.push_back(pml);
					auto boundary = Boundary::FindBoundary(pml, partition, partition->boundary_absorption_);
					boundaries_.push_back(boundary);
					info_.num_pml_partitions++;
					started = false;
				}
			}
			else if (partition->right_free_borders_[0][i])
			{
				start = i;
				end = i;
				started = true;
			}
		}

		// Find top PML.
		start = 0; end = 0; started = false;
		for (int i = 0; i < partition->width_; i++)
		{
			if (started)
			{
				end++;
				if (partition->top_free_borders_[0][i] && i != partition->width_ - 1)
				{
					continue;
				}
				else if (start != end)
				{
					if (i != partition->width_ - 1) end--;
					auto pml = std::make_shared<PmlPartition>(
						partition,
						PmlPartition::P_TOP,
						partition->x_start_ + start,
						partition->y_start_ - Simulation::n_pml_layers_,
						partition->z_start_,
						end - start + 1,
						Simulation::n_pml_layers_,
						partition->depth_);
					partitions_.push_back(pml);
					auto boundary = Boundary::FindBoundary(pml, partition, partition->boundary_absorption_);
					boundaries_.push_back(boundary);
					info_.num_pml_partitions++;
					started = false;
				}
			}
			else if (partition->top_free_borders_[0][i])
			{
				start = i;
				end = i;
				started = true;
			}
		}

		// Find bottom PML.
		start = 0; end = 0; started = false;
		for (int i = 0; i < partition->width_; i++)
		{
			if (started)
			{
				end++;
				if (partition->bottom_free_borders_[0][i] && i != partition->width_ - 1)
				{
					continue;
				}
				else if (start != end)
				{
					if (i != partition->width_ - 1) end--;
					auto pml = std::make_shared<PmlPartition>(
						partition,
						PmlPartition::P_BOTTOM,
						partition->x_start_ + start,
						partition->y_end_,
						partition->z_start_,
						end - start + 1,
						Simulation::n_pml_layers_,
						partition->depth_);
					partitions_.push_back(pml);
					auto boundary = Boundary::FindBoundary(pml, partition, partition->boundary_absorption_);
					boundaries_.push_back(boundary);
					info_.num_pml_partitions++;
					started = false;
				}
			}
			else if (partition->bottom_free_borders_[0][i])
			{
				start = i;
				end = i;
				started = true;
			}
		}

		// Add front PML.
		{
			auto pml = std::make_shared<PmlPartition>(
				partition,
				PmlPartition::P_FRONT,
				partition->x_start_,
				partition->y_start_,
				partition->z_start_ - Simulation::n_pml_layers_,
				partition->width_,
				partition->height_,
				Simulation::n_pml_layers_);
			partitions_.push_back(pml);
			std::shared_ptr<Boundary> cu_b(new Boundary(
				Boundary::Z_BOUNDARY,
				partition->boundary_absorption_,
				pml,
				partition,
				partition->x_start_,
				partition->x_end_,
				partition->y_start_,
				partition->y_end_,
				partition->z_start_ - 3,
				partition->z_start_ + 3));
			boundaries_.push_back(cu_b);
			info_.num_pml_partitions++;
		}


		// Add back PML.
		{
			auto pml = std::make_shared<PmlPartition>(
				partition,
				PmlPartition::P_BACK,
				partition->x_start_,
				partition->y_start_,
				partition->z_end_,
				partition->width_,
				partition->height_,
				Simulation::n_pml_layers_);
			partitions_.push_back(pml);
			std::shared_ptr<Boundary> cu_b(new Boundary(
				Boundary::Z_BOUNDARY,
				partition->boundary_absorption_,
				pml,
				partition,
				partition->x_start_,
				partition->x_end_,
				partition->y_start_,
				partition->y_end_,
				partition->z_end_ - 3,
				partition->z_end_ + 3));
			boundaries_.push_back(cu_b);
			info_.num_pml_partitions++;
		}
	}

	// Iterate over partitions_ and add elements that are not in dct_partitions_
	for (const auto& partition : partitions_) {
		// Check if the partition is not in dct_partitions_
		if (std::find(dct_partitions_.begin(), dct_partitions_.end(), partition) == dct_partitions_.end()) {
			pml_partitions_.push_back(partition); // Add to pml_partitions_ if not found in dct_partitions_
		}
	}

	/*------------- partitions includes pml partition --------------------------*/

	x_start_ = y_start_ = z_start_ = std::numeric_limits<int>::max();
	x_end_ = y_end_ = z_end_ = std::numeric_limits<int>::min();

	for (auto partition : partitions_)
	{
		x_start_ = std::min(x_start_, partition->x_start_);
		y_start_ = std::min(y_start_, partition->y_start_);
		z_start_ = std::min(z_start_, partition->z_start_);

		x_end_ = std::max(x_end_, partition->x_end_);
		y_end_ = std::max(y_end_, partition->y_end_);
		z_end_ = std::max(z_end_, partition->z_end_);
	}

	size_x_ = x_end_ - x_start_;
	size_y_ = y_end_ - y_start_;
	size_z_ = z_end_ - z_start_;

	// Each panel is square so all 3 fit neatly side-by-side
	panel_w_ = std::max({ size_x_, size_y_, size_z_ });
	panel_h_ = std::max({ size_x_, size_y_, size_z_ });

	pixels_.assign(panel_w_ * 3 * panel_h_, 0);
	cudaMalloc((void**)&d_pixels_, panel_w_ * 3 * panel_h_ * sizeof(uint32_t));
	sdl_fmt_ = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
	ready_ = true;
}

Simulation::~Simulation()
{
	if (sdl_fmt_) SDL_FreeFormat(sdl_fmt_);
	if (d_pixels_) cudaFree(d_pixels_);
}

int Simulation::Update()
{
	int time_step = time_step_++;
	//std::cout << "#" << std::setw(5) << time_step << " : ";
	//std::cout << std::to_string(sources_[0]->SampleValue(time_step)) << " ";

	// DEBUG: print pressure before first update
	if (time_step == 0 && !dct_partitions_.empty()) {
		double p = dct_partitions_[0]->get_pressure(5, 5, 3);
		std::cout << "Step 0 (Pre-Update) | p(5,5,3)=" << p << std::endl;
	}

	for (int i = 0; i < dct_partitions_.size(); i++)
	{
		// compute force
		dct_partitions_[i]->ComputeSourceForcingTerms((double)time_step);

		// update pressure and velocity
		dct_partitions_[i]->Update();

		// reset residue
		dct_partitions_[i]->reset_residues();
	}

	// DEBUG: print pressure values every 20 steps to diagnose stability
	if (time_step % 20 == 0 && time_step <= 400 && !dct_partitions_.empty()) {
		double p = dct_partitions_[0]->get_pressure(5, 5, 3);
		double src_val = sources_[0]->SampleValue(time_step);
		std::cout << "Step " << time_step << " | src=" << src_val
		          << " | p(5,5,3)=" << p << std::endl;
	}

	for (int i = 0; i < pml_partitions_.size(); i++)
	{
		// compute force
		pml_partitions_[i]->ComputeSourceForcingTerms((double)time_step);

		// update pressure and velocity
		pml_partitions_[i]->Update();

		// reset residue
		pml_partitions_[i]->reset_residues();
	}

	// Wait for all partition internal updates to finish before boundary calculations
	cudaDeviceSynchronize();

	// compute residue
	for (auto p : partitions_) p->reset_residues();

	if (time_step == 0) {
		std::cout << "Runtime: updating " << boundaries_.size() << " boundaries." << std::endl;
	}

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < boundaries_.size(); i++) {
		boundaries_[i]->ComputeResidues();
	}
	// Post-merge depends on residues
	cudaDeviceSynchronize();

	// post-merge
	for (int i = 0; i < partitions_.size(); i++)
	{
		partitions_[i]->PostMerge();
	}
	
	// Clear forces for next step (can be async)
	for (auto p : partitions_) p->reset_forces();

	cudaDeviceSynchronize();
	//std::cout << std::endl;

	// Visualization: render XY / XZ / YZ planes side-by-side every viz_skip_ steps
	if (time_step % viz_skip_ == 0)
	{
		float v_coef = 0.1f;

		// Anchor slices through the first source position
		int src_x = sources_[0]->x() - x_start_;
		int src_y = sources_[0]->y() - y_start_;
		int src_z = sources_[0]->z() - z_start_;

		cudaMemset(d_pixels_, 0, panel_w_ * 3 * panel_h_ * sizeof(uint32_t));

		for (auto partition : partitions_)
		{
			int x_off = partition->x_start_ - x_start_;
			int y_off = partition->y_start_ - y_start_;
			int z_off = partition->z_start_ - z_start_;

			// Panel 0: XY plane (constant Z = src_z)
			// x→screen-x, y→screen-y, offset=(0, 0)
			partition->RenderToBuffer(d_pixels_,
				0, src_z + z_start_,       // plane_type=0(XY), global z coord
				panel_w_ * 3, panel_h_,    // full buffer width, panel height
				x_off, y_off, v_coef);

			// Panel 1: XZ plane (constant Y = src_y)
			// x→screen-x, z→screen-y, offset=(panel_w_, 0)
			partition->RenderToBuffer(d_pixels_,
				2, src_y + y_start_,       // plane_type=2(XZ), global y coord
				panel_w_ * 3, panel_h_,
				panel_w_ + x_off, z_off, v_coef);

			// Panel 2: YZ plane (constant X = src_x)
			// y→screen-x, z→screen-y, offset=(2*panel_w_, 0)
			partition->RenderToBuffer(d_pixels_,
				1, src_x + x_start_,       // plane_type=1(YZ), global x coord
				panel_w_ * 3, panel_h_,
				2 * panel_w_ + y_off, z_off, v_coef);
		}

		// Single PCIe transfer for the entire 3-panel frame
		cudaMemcpy(pixels_.data(), d_pixels_,
			panel_w_ * 3 * panel_h_ * sizeof(uint32_t), cudaMemcpyDeviceToHost);
	}

	return time_step;
}

void Simulation::Info()
{
	std::cout << "# Simulation Info. #########################################" << std::endl;
	std::cout << "Simulation: "
		<< std::to_string(x_start_) << "," << std::to_string(y_start_) << "," << std::to_string(z_start_) << "->"
		<< std::to_string(x_end_) << "," + std::to_string(y_end_) << "," + std::to_string(z_end_) << std::endl;
	std::cout << "Size: " << "<" << size_x_ << "," << size_y_ << "," << size_z_ << ">" << std::endl;
	std::cout << "dh = " << std::to_string(Simulation::dh_)
		<< "(m), dt = " << std::to_string(Simulation::dt_)
		<< "(s), c0 = " << std::to_string(Simulation::c0_)
		<< "(m/s), alpha1 = " << std::to_string(Simulation::air_absorption_alpha1_)
		<< "(1/s), alpha2 = " << std::to_string(Simulation::air_absorption_alpha2_) << "(1/s)"
		<< std::endl;
	std::cout << "Number of dct_partitions: " << info_.num_dct_partitions << std::endl;
	std::cout << "Number of pml_partitions: " << info_.num_pml_partitions << std::endl;
	std::cout << "Number of boundaries: " << info_.num_boundaries << std::endl;
	std::cout << "Number of sources: " << info_.num_sources << std::endl;

	std::cout << "############################################################" << std::endl;
	for (auto p : partitions_)
	{
		if (p->info_.type == "DCT")
			p->Info();
	}
	std::cout << "------------------------------------------------------------" << std::endl;
	//for (auto b : boundaries_)
	//{
	//	b->Info();
	//}
	std::cout << "------------------------------------------------------------" << std::endl;
	for (auto s : sources_)
	{
		std::cout << "Source " << s->id_ << ": " << s->x() << "," << s->y() << "," << s->z() << std::endl;
	}
	std::cout << "------------------------------------------------------------" << std::endl;
}

