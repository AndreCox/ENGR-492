#include "spring_system.h"

SpringSystem::SpringSystem(std::vector<Node> &n, std::vector<Spring> &s)
    : nodes(n), springs(s), max_stress(0.0f), min_stress(0.0f)
{
    // Resize displacement and forces vectors
    int total_dof = static_cast<int>(nodes.size()) * 2; // 2 DOF per node (x and y)
    displacement = Eigen::VectorXd::Zero(total_dof);
    forces = Eigen::VectorXd::Zero(total_dof);

    // step 1 - compute stiffness matrices for each spring
    for (auto &spring : springs)
    {
        spring.compute_stiffness(nodes);
    }

    // step 2 - assemble global stiffness matrix
    assemble_global_stiffness();
}

// MPC (Multi-Point Constraint) for slider nodes
// This generates a constraint equation: a_x * u + a_y * v = 0
// which means displacement perpendicular to the slider direction is zero
void SpringSystem::generate_constraint_row(Eigen::MatrixXd &C, int row_index, int node_id, const std::vector<int> &free_dof_indices)
{
    const Node &node = nodes[node_id];
    if (node.constraint_type != Slider)
        return; // MPC only for slider nodes

    // Slider angle in radians
    double theta = node.constraint_angle * static_cast<double>(M_PI) / 180.0f;
    double normal_angle = theta + static_cast<double>(M_PI) / 2.0f;

    double a_x = std::cos(normal_angle);
    double a_y = std::sin(normal_angle);

    int dof_x = node.id * 2;
    int dof_y = node.id * 2 + 1;

    // this is essentially calculating a_x * u_dof + a_y * v_dof = 0

    // Find positions of these DOFs in the reduced system and set entries
    for (int j = 0; j < static_cast<int>(free_dof_indices.size()); ++j)
    {
        if (free_dof_indices[j] == dof_x)
        {
            C(row_index, j) = a_x;
        }
        else if (free_dof_indices[j] == dof_y)
        {
            C(row_index, j) = a_y;
        }
    }
}

// step 3 - attempt to solve the system
int SpringSystem::solve_system()
{
    int num_nodes = static_cast<int>(nodes.size());
    int total_dof = num_nodes * 2; // 2 DOF per node (x and y)

    // Identify free DOFs (not fixed nodes)
    std::vector<int> free_dof_indices;
    free_dof_indices.reserve(total_dof);

    for (int i = 0; i < num_nodes; ++i)
    {
        if (nodes[i].constraint_type != Fixed)
        {
            free_dof_indices.push_back(i * 2);     // x DOF
            free_dof_indices.push_back(i * 2 + 1); // y DOF
        }
    }

    int num_free_dofs = static_cast<int>(free_dof_indices.size());

    if (num_free_dofs == 0)
    {
        std::cout << "No free DOFs to solve!" << std::endl;
        return -1;
    }

    // Create reduced stiffness matrix and force vector
    Eigen::MatrixXd K_r = Eigen::MatrixXd::Zero(num_free_dofs, num_free_dofs);
    Eigen::VectorXd F_r = Eigen::VectorXd::Zero(num_free_dofs);

    for (int i = 0; i < num_free_dofs; ++i)
    {
        for (int j = 0; j < num_free_dofs; ++j)
        {
            K_r(i, j) = global_k_matrix(free_dof_indices[i], free_dof_indices[j]);
        }
        F_r(i) = forces(free_dof_indices[i]);
    }

    std::cout << "Reduced Stiffness Matrix (" << num_free_dofs << "x" << num_free_dofs << "):\n"
              << K_r << std::endl;
    std::cout << "Reduced Force Vector:\n"
              << F_r << std::endl;

    // Identify slider nodes and create constraint matrix
    std::vector<int> slider_nodes;
    for (int i = 0; i < num_nodes; ++i)
    {
        if (nodes[i].constraint_type == Slider)
            slider_nodes.push_back(i);
    }

    int num_constraints = static_cast<int>(slider_nodes.size());

    if (num_constraints == 0)
    {
        // No constraints - solve directly
        Eigen::VectorXd u_r = K_r.fullPivLu().solve(F_r);

        displacement = Eigen::VectorXd::Zero(total_dof);
        for (int i = 0; i < num_free_dofs; ++i)
        {
            displacement(free_dof_indices[i]) = u_r(i);
        }
    }
    else
    {
        // Build constraint matrix in reduced coordinates using generate_constraint_row
        Eigen::MatrixXd C_r = Eigen::MatrixXd::Zero(num_constraints, num_free_dofs);

        for (int ci = 0; ci < num_constraints; ++ci)
        {
            int node_id = slider_nodes[ci];
            generate_constraint_row(C_r, ci, node_id, free_dof_indices);
        }

        std::cout << "Constraint Matrix C_r (" << num_constraints << "x" << num_free_dofs << "):\n"
                  << C_r << std::endl;

        // Scale the constraints to match the stiffness matrix magnitude for better conditioning
        double k_scale = K_r.norm();
        double constraint_scale = (k_scale > 0.0) ? k_scale : 1.0;
        Eigen::MatrixXd C_r_scaled = C_r * constraint_scale;

        // Create augmented saddle point system
        int augmented_size = num_free_dofs + num_constraints;
        Eigen::MatrixXd saddle_matrix = Eigen::MatrixXd::Zero(augmented_size, augmented_size);
        Eigen::VectorXd saddle_rhs = Eigen::VectorXd::Zero(augmented_size);

        // Fill in K_r block
        saddle_matrix.block(0, 0, num_free_dofs, num_free_dofs) = K_r;

        // Fill in scaled C_r^T block (upper right)
        saddle_matrix.block(0, num_free_dofs, num_free_dofs, num_constraints) = C_r_scaled.transpose();

        // Fill in scaled C_r block (lower left)
        saddle_matrix.block(num_free_dofs, 0, num_constraints, num_free_dofs) = C_r_scaled;

        // Fill in RHS
        saddle_rhs.head(num_free_dofs) = F_r;

        // print out the saddle matrix and rhs for debugging
        std::cout << "Saddle Point Matrix (" << augmented_size << "x" << augmented_size << "):\n"
                  << saddle_matrix << std::endl;
        std::cout << "Saddle Point RHS:\n"
                  << saddle_rhs << std::endl;

        // Solve the saddle point system
        Eigen::VectorXd full_solution = saddle_matrix.fullPivLu().solve(saddle_rhs);

        if (full_solution.size() != augmented_size)
        {
            std::cerr << "Saddle point solver failed or returned unexpected size." << std::endl;
            return -2;
        }

        // Extract reduced displacements and (scaled) Lagrange multipliers
        Eigen::VectorXd u_r = full_solution.head(num_free_dofs);
        Eigen::VectorXd lagrange_multipliers_scaled = full_solution.tail(num_constraints);

        // Unscale Lagrange multipliers back
        Eigen::VectorXd lagrange_multipliers = lagrange_multipliers_scaled * constraint_scale;

        // Map back to global displacement vector
        displacement = Eigen::VectorXd::Zero(total_dof);
        for (int i = 0; i < num_free_dofs; ++i)
        {
            displacement(free_dof_indices[i]) = u_r(i);
        }
    }

    // Print solution summary
    std::cout << "\n=== SOLUTION ===\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        double u = displacement(i * 2);
        double v = displacement(i * 2 + 1);
        double total_disp = std::sqrt(u * u + v * v);
        std::cout << "  Node " << i << ": u=" << u << " m, v=" << v << " m (total=" << total_disp << " m)\n";

        if (nodes[i].constraint_type == Slider)
        {
            double theta = nodes[i].constraint_angle * static_cast<double>(M_PI) / 180.0f;
            double along_slider = u * std::cos(theta) + v * std::sin(theta);
            double perp_slider = -u * std::sin(theta) + v * std::cos(theta);
            std::cout << "    Movement along slider (" << nodes[i].constraint_angle << "°): " << along_slider << " m\n";
            std::cout << "    Movement perpendicular to slider: " << perp_slider << " m (should be ~0)\n";
        }
    }

    // Calculate reaction forces at fixed and slider nodes
    Eigen::VectorXd reactions = global_k_matrix * displacement - forces;
    std::cout << "\nReaction Forces (N):\n";

    double total_reaction_x = 0.0;
    double total_reaction_y = 0.0;

    for (int i = 0; i < num_nodes; ++i)
    {
        if (nodes[i].constraint_type != Free)
        {
            std::string constraint_str = (nodes[i].constraint_type == Fixed) ? "Fixed" : "Slider";
            double rx = reactions(i * 2);
            double ry = reactions(i * 2 + 1);
            std::cout << "  Node " << i << " (" << constraint_str << "): Fx=" << rx << " N, Fy=" << ry << " N\n";
            total_reaction_x += rx;
            total_reaction_y += ry;
        }
    }

    // Equilibrium check
    double total_applied_x = 0.0;
    double total_applied_y = 0.0;
    for (int i = 0; i < num_nodes; ++i)
    {
        total_applied_x += forces(i * 2);
        total_applied_y += forces(i * 2 + 1);
    }

    std::cout << "\nEquilibrium Check:\n";
    std::cout << "  Applied Fx = " << total_applied_x << " N\n";
    std::cout << "  Applied Fy = " << total_applied_y << " N\n";
    std::cout << "  Reaction Fx = " << total_reaction_x << " N\n";
    std::cout << "  Reaction Fy = " << total_reaction_y << " N\n";
    std::cout << "  Balance (Fx): " << (total_applied_x + total_reaction_x) << " N (should be ~0)\n";
    std::cout << "  Balance (Fy): " << (total_applied_y + total_reaction_y) << " N (should be ~0)\n";

    // Calculate internal forces and stresses in springs
    std::cout << "\nSpring Internal Forces (N, tension positive):\n";
    max_stress = -1e10f;
    min_stress = 1e10f;

    for (auto &spring : springs)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        Eigen::Vector4d element_disp;
        element_disp << displacement(n1 * 2), displacement(n1 * 2 + 1),
            displacement(n2 * 2), displacement(n2 * 2 + 1);

        Eigen::Vector4d element_forces = spring.k_matrix * element_disp;

        double dx = static_cast<double>(nodes[n2].position[0] - nodes[n1].position[0]);
        double dy = static_cast<double>(nodes[n2].position[1] - nodes[n1].position[1]);
        double length = std::sqrt(dx * dx + dy * dy);
        double c = dx / length;
        double s = dy / length;

        // Axial force (negative of force at first node)
        double axial_force = -(c * element_forces(0) + s * element_forces(1));

        // Stress = Force / Area (convert to MPa if area in mm^2 — keep consistent with your units)
        spring.stress = static_cast<double>(axial_force / spring.A / 1e6); // MPa

        max_stress = std::max(max_stress, spring.stress);
        min_stress = std::min(min_stress, spring.stress);

        std::cout << "  Spring " << spring.id << " (nodes " << n1 << "-" << n2 << "): " << axial_force << " N\n";
    }

    std::cout << "\nSpring Axial Stresses (MPa):\n";
    for (const auto &spring : springs)
    {
        std::cout << "  Spring " << spring.id << " (nodes " << spring.nodes[0] << "-" << spring.nodes[1]
                  << "): " << spring.stress << " MPa\n";
    }

    std::cout << "\nStress Range: " << min_stress << " to " << max_stress << " MPa\n";

    return 0;
}

void SpringSystem::assemble_global_stiffness()
{
    int num_nodes = static_cast<int>(nodes.size());
    int total_dof = num_nodes * 2;
    global_k_matrix = Eigen::MatrixXd::Zero(total_dof, total_dof);

    // Assemble contributions from each spring
    for (const auto &spring : springs)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        int dof1_x = n1 * 2;
        int dof1_y = n1 * 2 + 1;
        int dof2_x = n2 * 2;
        int dof2_y = n2 * 2 + 1;

        int dofs[4] = {dof1_x, dof1_y, dof2_x, dof2_y};

        // Add element stiffness matrix to global matrix
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                global_k_matrix(dofs[i], dofs[j]) += spring.k_matrix(i, j);
            }
        }
    }

    std::cout << "Global Stiffness Matrix (" << total_dof << "x" << total_dof << "):\n"
              << global_k_matrix << std::endl;
}
