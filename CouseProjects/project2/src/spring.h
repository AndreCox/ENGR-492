#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "node.h"
#include "material_profile.h"

class Spring
{
public:
    Spring(int n1, int n2, BeamProperties &properties)
        : properties(properties)
    {
        nodes[0] = n1;
        nodes[1] = n2;
        stress = 0.0f;
    }

    int nodes[2]; // node indices
    double k;     // stiffness (EA/L)
    float stress; // axial stress in MPa

    BeamProperties properties; // material and cross-sectional properties

    Eigen::Matrix4d k_matrix; // 4x4 stiffness matrix in global coordinates
                              // [u1, v1, u2, v2]

    // step 1 - compute the stiffness matrix for the spring
    void compute_stiffness(const std::vector<Node> &node_list);
};