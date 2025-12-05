#pragma once

#include <vector>
#include <Eigen/Eigen>
#include "node.h"
#include "spring.h"
#include "material_profile.h"
#include <iostream>
#include <cmath>

class SpringSystem
{
public:
    SpringSystem(std::vector<Node> &n, std::vector<Spring> &s);

    void generate_constraint_row(Eigen::MatrixXd &C, int row_index, int node_id, const std::vector<int> &free_dof_indices);
    int solve_system();
    void assemble_global_stiffness();

    std::vector<Node> nodes;
    std::vector<Spring> springs;
    std::vector<BeamProperties> spring_properties;
    Eigen::MatrixXd global_k_matrix;
    Eigen::VectorXd forces;       // forces in x and y directions [F1x, F1y, F2x, F2y, ...]
    Eigen::VectorXd displacement; // displacements in x and y [u1, v1, u2, v2, ...]
    int total_dof;
    float max_stress;
    float min_stress;
};
