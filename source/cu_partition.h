#pragma once
#include "partition.h"
#include <vector>
#include <string>

class CuBoundary;

// CuPartition extends Partition so it can be stored in the existing
// vector<shared_ptr<Partition>> throughout simulation.cpp and partition.cpp.
// The heavy data lives on VRAM (d_pressure_, d_velocity_, etc.) while 
// the CPU-side virtual accessor methods perform device-to-host copies on demand.
class CuPartition : public Partition
{
public:
    CuPartition(int xs, int ys, int zs, int w, int h, int d);
    virtual ~CuPartition();

    // GPU-resident data pointers (device memory)
    double* d_pressure_;
    double* d_velocity_;
    double* d_force_;
    double* d_residue_;

    bool should_render_ = true;

    virtual void Update() = 0;

    // --- Partition virtual interface implementation ---
    // These perform host-device copies so that the existing CPU simulation
    // wiring (boundaries, sources, rendering) can access data transparently.
    virtual double* get_pressure_field() override;

    virtual double get_pressure(int x, int y, int z) override;
    virtual void set_pressure(int x, int y, int z, double v) override;
    virtual void add_to_pressure(int x, int y, int z, double v) override;

    virtual double get_velocity(int x, int y, int z) override;
    virtual void set_velocity(int x, int y, int z, double v) override;
    virtual void add_to_velocity(int x, int y, int z, double v) override;

    virtual double get_residue(int x, int y, int z) override;
    virtual void set_residue(int x, int y, int z, double v) override;
    virtual void add_to_residue(int x, int y, int z, double v) override;

    virtual double get_force(int x, int y, int z) override;
    virtual void set_force(int x, int y, int z, double v) override;

    virtual void reset_forces() override;
    virtual void reset_residues() override;

    virtual void PostMerge() override;
    virtual void ComputeSourceForcingTerms(double t) override;

    // --- GPU-side Visualization and Source Forcing ---
    // d_pixels: buffer to render into (RGBA8888)
    // plane_type: 0=XY, 1=YZ, 2=XZ
    // coord: the z, x, or y coordinate respectively
    void RenderToBuffer(uint32_t* d_pixels, int plane_type, int coord, int screen_width, int screen_height, int x_offset, int y_offset, float v_coef);

    // Batched sources
    int num_sources_ = 0;
    int* d_source_indices_ = nullptr; // flat indices: z * H * W + y * W + x
    double* d_source_values_ = nullptr; // current time-step values

    virtual std::vector<double> get_xy_plane(int z) override;
    virtual std::vector<double> get_yz_plane(int x) override;
    virtual std::vector<double> get_xz_plane(int y) override;

    // Boundary handling (CUDA version)
    void AddCuBoundary(CuBoundary* b);
    std::vector<CuBoundary*> cu_boundaries_;
};
