#pragma once

#include <vector>
#include <Eigen/Eigen>
#include "node.h"
#include "beam.h"
#include "beam_props.h"
#include <iostream>
#include <cmath>

enum UnitSystem
{
    Imperial,
    Metric
};

class FEMSystem
{
public:
    FEMSystem(std::vector<Node> &n, std::vector<Beam> &s, std::vector<MaterialProfile> &materials, std::vector<BeamProfile> &beam_profiles);

    UnitSystem unit_system = Imperial;

    void generate_constraint_row(Eigen::MatrixXd &C, int row_index, int node_id, const std::vector<int> &free_dof_indices);
    int solve_system();
    void assemble_global_stiffness();

    std::vector<Node> nodes;
    std::vector<Beam> beams;
    std::vector<MaterialProfile> materials_list;
    std::vector<BeamProfile> beam_profiles_list;
    Eigen::MatrixXd global_k_matrix;
    Eigen::VectorXd forces;       // forces in x and y directions [F1x, F1y, F2x, F2y, ...]
    Eigen::VectorXd displacement; // displacements in x and y [u1, v1, u2, v2, ...]
    int total_dof;
    float max_stress;
    float min_stress;
};
