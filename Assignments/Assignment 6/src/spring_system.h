#pragma once

#include <vector>
#include <Eigen/Eigen>
#include "node.h"
#include "spring.h"
#include <iostream>
#include <cmath>

class SpringSystem
{
public:
    // Only declare the constructor here (no inline body)
    SpringSystem(const std::vector<Node> &n, const std::vector<Spring> &s);

    void generate_constraint_row(Eigen::MatrixXd &C, int row_index, int node_id, const std::vector<int> &free_dof_indices);
    int solve_system();
    void assemble_global_stiffness();

    std::vector<Node> nodes;
    std::vector<Spring> springs;
    Eigen::MatrixXd global_k_matrix;
    Eigen::VectorXd forces;       // forces in x and y directions [F1x, F1y, F2x, F2y, ...]
    Eigen::VectorXd displacement; // displacements in x and y [u1, v1, u2, v2, ...]
    float max_stress;
    float min_stress;
};
