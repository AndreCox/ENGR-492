#include "spring.h"
// step 1 - compute the stiffness matrix for the spring
void Spring::compute_stiffness(const std::vector<Node> &node_list)
{
    int n1 = nodes[0];
    int n2 = nodes[1];

    float dx = node_list[n2].position[0] - node_list[n1].position[0];
    float dy = node_list[n2].position[1] - node_list[n1].position[1];
    float length = std::sqrt(dx * dx + dy * dy);

    // EA/L stiffness
    k = (A * E) / length;

    // Direction cosines
    float c = dx / length; // cos(theta)
    float s = dy / length; // sin(theta)

    // Local stiffness matrix for a 2D truss element
    float c2 = c * c;
    float s2 = s * s;
    float cs = c * s;

    k_matrix << c2, cs, -c2, -cs,
        cs, s2, -cs, -s2,
        -c2, -cs, c2, cs,
        -cs, -s2, cs, s2;

    k_matrix *= k;
}
