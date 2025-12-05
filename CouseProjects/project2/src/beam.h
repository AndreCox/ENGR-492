#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "node.h"
#include "beam_props.h"

class Beam
{
public:
    int nodes[2];
    double k;
    double axial_force;
    double max_moment;
    float stress;

    int material_idx;
    int shape_idx;

    Eigen::MatrixXd k_matrix;

    Beam() : nodes{-1, -1}, k(0.0), stress(0.0f), material_idx(-1), shape_idx(-1)
    {
        k_matrix.setZero();
    }

    Beam(int n1, int n2, int mat, int shp)
        : nodes{n1, n2}, k(0.0), stress(0.0f),
          material_idx(mat), shape_idx(shp)
    {
        k_matrix.setZero();
    }

    void compute_stiffness(const std::vector<Node> &node_list,
                           const std::vector<MaterialProfile> &materials,
                           const std::vector<BeamProfile> &shapes);
};