#pragma once
#include "partition.h"

class PmlPartition : public Partition
{
public:
    enum PmlType {
        P_LEFT, P_RIGHT, P_TOP, P_BOTTOM, P_FRONT, P_BACK
    } type_;

    // Auxiliary PML fields (device memory)
    double* d_psi_;
    double* d_phi_x_;
    double* d_phi_y_;
    double* d_phi_z_;
    double* d_zetax_;
    double* d_zetay_;
    double* d_zetaz_;

    // To avoid creating new buffers every step, we use temp pointers for swaps
    double* d_p_new_;
    double* d_v_new_;
    double* d_psi_new_;
    double* d_phi_x_new_;
    double* d_phi_y_new_;
    double* d_phi_z_new_;

    PmlPartition(std::shared_ptr<Partition> neighbor_part, PmlType type, int xs, int ys, int zs, int w, int h, int d);
    virtual ~PmlPartition();

    virtual void Update() override;
    
    // We reuse CuPartition's data accessors but we might need to sync more buffers if debugging
};
