#include "simulation.h"
#include "partition.h"
#include "pml_partition.h"
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
			boundaries_.push_back(std::make_shared<Boundary>(
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
			boundaries_.push_back(std::make_shared<Boundary>(
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



	pixels_.assign(size_x_*size_y_, 0);
	ready_ = true;
}

Simulation::~Simulation()
{
}

int Simulation::Update()
{
	int time_step = time_step_++;
	//std::cout << "#" << std::setw(5) << time_step << " : ";
	//std::cout << std::to_string(sources_[0]->SampleValue(time_step)) << " ";

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < dct_partitions_.size(); i++)
	{
		// reset force
		//dct_partitions_[i]->reset_forces();

		// compute force
		dct_partitions_[i]->ComputeSourceForcingTerms(time_step);
		//std::cout << "impose force partition " << partition->info_.id << " ";

		// update pressure and velocity
		dct_partitions_[i]->Update();
		//std::cout << "update pressure partition " << partition->info_.id << " ";

		// reset residue
		dct_partitions_[i]->reset_residues();
	}

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < pml_partitions_.size(); i++)
	{
		// reset force
		//pml_partitions_[i]->reset_forces();

		// compute force
		pml_partitions_[i]->ComputeSourceForcingTerms(time_step);
		//std::cout << "impose force partition " << partition->info_.id << " ";

		// update pressure and velocity
		pml_partitions_[i]->Update();
		//std::cout << "update pressure partition " << partition->info_.id << " ";

		// reset residue
		pml_partitions_[i]->reset_residues();
	}

	// compute residue
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < boundaries_.size(); i++) {
		boundaries_[i]->ComputeResidues();
	}

	//std::cout << std::endl;

	// post-merge
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < partitions_.size(); i++)
	{
		partitions_[i]->PostMerge();
	}
	//std::cout << std::endl;

	// Visualization
	if (time_step % 1 == 0)
	{
		SDL_PixelFormat* fmt = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
		double v_coef = 1.0;
		bool render_pml = true;
		if (look_from_ == 0)	//xy
		{
			int pixels_z = sources_[0]->z();
			for (auto partition : partitions_)
			{
				if (!render_pml)
				{
					if (!partition->should_render_) continue;
				}
				//if (partition->is_z_pml_) continue;
				if (partition->z_start_ > pixels_z || partition->z_end_ < pixels_z)
				{
					continue;
				}
				int x_offset = partition->x_start_ - x_start_;
				int y_offset = partition->y_start_ - y_start_;
				std::vector<double> partition_xy;
				partition_xy = partition->get_xy_plane(pixels_z);

				int height = partition->height_;
				int width = partition->width_;

#pragma omp parallel for schedule(dynamic)
				for (int idx = 0; idx < height * width; idx++) {
					// Compute the indices i and j from the flattened index idx
					int i = idx / width;
					int j = idx % width;

					double pressure = partition_xy[i * width + j];
					double norm = 0.5 * std::max(-1.0, std::min(1.0, pressure * v_coef)) + 0.5;
					int r, g, b;

					if (norm >= 0.5)
					{
						r = static_cast<int>(255 - round(255.0 * 2.0 * (norm - 0.5)));
						g = static_cast<int>(255 - round(255.0 * 2.0 * (norm - 0.5)));
						b = 255;
					}
					else {
						r = 255;
						g = static_cast<int>(255 - round(255.0 * (1.0 - 2.0 * norm)));
						b = static_cast<int>(255 - round(255.0 * (1.0 - 2.0 * norm)));
					}

					if (partition->should_render_)
					{
						pixels_[(y_offset + i) * size_x_ + (x_offset + j)] = SDL_MapRGBA(fmt, 255, r, g, b);
					}
					else
					{
						pixels_[(y_offset + i) * size_x_ + (x_offset + j)] = SDL_MapRGBA(fmt, 255, 0.5 * r, 0.5 * g, 0.5 * b);
					}
				}
			}
		}
		else if (look_from_ == 1)	//yz
		{
			int pixels_x = sources_[0]->x();
			for (auto partition : partitions_)
			{
				if (!render_pml)
				{
					if (!partition->should_render_) continue;
				}
				//if (partition->is_x_pml_) continue;
				if (partition->x_start_ > pixels_x || partition->x_end_ < pixels_x)
				{
					continue;
				}
				int y_offset = partition->y_start_ - y_start_;
				int z_offset = partition->z_start_ - z_start_;
				std::vector<double> partition_yz;
				partition_yz = partition->get_yz_plane(pixels_x);

				int depth = partition->depth_;
				int height = partition->height_;

#pragma omp parallel for
				for (int idx = 0; idx < depth * height; idx++)
				{
					// Compute the original indices i and j from the flattened index
					int i = idx / height;
					int j = idx % height;

					double pressure = partition_yz[idx];
					double norm = 0.5 * std::max(-1.0, std::min(1.0, pressure * v_coef)) + 0.5;
					int r, g, b;

					if (norm >= 0.5)
					{
						r = static_cast<int>(255 - round(255.0 * 2.0 * (norm - 0.5)));
						g = static_cast<int>(255 - round(255.0 * 2.0 * (norm - 0.5)));
						b = 255;
					}
					else
					{
						r = 255;
						g = static_cast<int>(255 - round(255.0 * (1.0 - 2.0 * norm)));
						b = static_cast<int>(255 - round(255.0 * (1.0 - 2.0 * norm)));
					}

					if (partition->should_render_)
					{
						pixels_[(z_offset + i) * size_y_ + (y_offset + j)] = SDL_MapRGBA(fmt, 255, r, g, b);
					}
					else
					{
						pixels_[(z_offset + i) * size_y_ + (y_offset + j)] = SDL_MapRGBA(fmt, 255, 0.5 * r, 0.5 * g, 0.5 * b);
					}
				}
			}
		}
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

