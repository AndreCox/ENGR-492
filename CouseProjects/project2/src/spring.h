#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "node.h"

class Spring
{
public:
    Spring(int i, int n1, int n2, double A, double E) : id(i), A(A), E(E)
    {
        nodes[0] = n1;
        nodes[1] = n2;
        stress = 0.0f;
    }

    int id;
    int nodes[2]; // node indices
    double k;     // stiffness (EA/L)
    float stress; // axial stress in MPa

    double A;
    double E;

    Eigen::Matrix4d k_matrix; // 4x4 stiffness matrix in global coordinates
                              // [u1, v1, u2, v2]

    // step 1 - compute the stiffness matrix for the spring
    void compute_stiffness(const std::vector<Node> &node_list);
};