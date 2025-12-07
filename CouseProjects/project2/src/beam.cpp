#include "beam.h"
// step 1 - compute the stiffness matrix for the spring
void Beam::compute_stiffness(const std::vector<Node> &node_list,
                             const std::vector<MaterialProfile> &materials,
                             const std::vector<BeamProfile> &shapes)
{
    int n1 = nodes[0];
    int n2 = nodes[1];

    const MaterialProfile &material = materials[material_idx];
    const BeamProfile &shape = shapes[shape_idx];

    k_matrix.resize(6, 6);

    double E = material.youngs_modulus;
    double A = shape.area;
    double I = shape.moment_of_inertia;

    if (is_truss)
    {
        I = 0.0; // For truss elements, set moment of inertia to zero
    }

    double dx = node_list[n2].position[0] - node_list[n1].position[0];
    double dy = node_list[n2].position[1] - node_list[n1].position[1];
    double L = std::sqrt(dx * dx + dy * dy);

    if (L < 1e-9)
    {
        k_matrix.setZero(); // Handle zero-length element
        return;
    }

    // EA/L stiffness
    k = (shape.area * material.youngs_modulus) / L;

    // Direction cosines
    double c = dx / L; // cos(theta)
    double s = dy / L; // sin(theta)

    // Calculate terms for stiffness matrix
    double EA_over_L = (E * A) / L;
    double EI_over_L = (E * I) / L;
    double EI_over_L2 = EI_over_L / L;
    double EI_over_L3 = EI_over_L2 / L;

    Eigen::MatrixXd k_prime = Eigen::MatrixXd::Zero(6, 6);

    // Axial (Row/Col 0 and 3)
    k_prime(0, 0) = EA_over_L;
    k_prime(0, 3) = -EA_over_L;
    k_prime(3, 0) = -EA_over_L;
    k_prime(3, 3) = EA_over_L;

    // Bending/Shear (Rows/Cols 1, 2, 4, 5)
    // Row 1 (V1)
    k_prime(1, 1) = 12.0 * EI_over_L3;
    k_prime(1, 2) = 6.0 * EI_over_L2;
    k_prime(1, 4) = -12.0 * EI_over_L3;
    k_prime(1, 5) = 6.0 * EI_over_L2;

    // Row 2 (M1)
    k_prime(2, 1) = 6.0 * EI_over_L2;
    k_prime(2, 2) = 4.0 * EI_over_L;
    k_prime(2, 4) = -6.0 * EI_over_L2;
    k_prime(2, 5) = 2.0 * EI_over_L;

    // Row 4 (V2)
    k_prime(4, 1) = -12.0 * EI_over_L3;
    k_prime(4, 2) = -6.0 * EI_over_L2;
    k_prime(4, 4) = 12.0 * EI_over_L3;
    k_prime(4, 5) = -6.0 * EI_over_L2;

    // Row 5 (M2)
    k_prime(5, 1) = 6.0 * EI_over_L2;
    k_prime(5, 2) = 2.0 * EI_over_L;
    k_prime(5, 4) = -6.0 * EI_over_L2;
    k_prime(5, 5) = 4.0 * EI_over_L;

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);

    T(0, 0) = c;
    T(0, 1) = s;
    T(1, 0) = -s;
    T(1, 1) = c;

    T(2, 2) = 1.0; // rotation about z-axis remains unchanged

    T(3, 3) = c;
    T(3, 4) = s;
    T(4, 3) = -s;
    T(4, 4) = c;

    T(5, 5) = 1.0; // rotation about z-axis remains unchanged

    k_matrix = T.transpose() * k_prime * T;
}
