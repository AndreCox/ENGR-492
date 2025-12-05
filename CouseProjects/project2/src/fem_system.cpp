#include "fem_system.h"

FEMSystem::FEMSystem(std::vector<Node> &n, std::vector<Beam> &s, std::vector<MaterialProfile> &materials, std::vector<BeamProfile> &beam_profiles)

    : nodes(n), beams(s), materials_list(materials), beam_profiles_list(beam_profiles), max_stress(0.0f), min_stress(0.0f)
{
    // Resize displacement and forces vectors
    total_dof = static_cast<int>(nodes.size()) * 3; // 3 DOF per node (x,y and theta)
    displacement = Eigen::VectorXd::Zero(total_dof);
    forces = Eigen::VectorXd::Zero(total_dof);
}

// MPC (Multi-Point Constraint) for slider nodes
// This generates a constraint equation: a_x * u + a_y * v = 0
// which means displacement perpendicular to the slider direction is zero
void FEMSystem::generate_constraint_row(Eigen::MatrixXd &C, int row_index, int index, const std::vector<int> &free_dof_indices)
{
    const Node &node = nodes[index];
    if (node.constraint_type != Slider)
        return; // MPC only for slider nodes

    // Slider angle in radians
    double theta = node.constraint_angle * static_cast<double>(M_PI) / 180.0f;
    double normal_angle = theta + static_cast<double>(M_PI) / 2.0f;

    double a_x = std::cos(normal_angle);
    double a_y = std::sin(normal_angle);

    int dof_x = index * 3;
    int dof_y = index * 3 + 1;
    // theta DOF is not used as rotation about slider is free

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

//  attempt to solve the system
int FEMSystem::solve_system()
{
    int num_nodes = static_cast<int>(nodes.size());
    if (num_nodes * 3 != total_dof) // check if total_dof needs updating
    {

        total_dof = static_cast<int>(nodes.size()) * 3; // 2 DOF per node (x and y)
        // shift the existing vectors to match new size
        displacement.conservativeResize(total_dof);
        forces.conservativeResize(total_dof);
    }

    // step 1 - compute stiffness matrices for each spring
    for (auto &beam : beams)
    {
        beam.compute_stiffness(nodes, materials_list, beam_profiles_list);
    }

    // step 2 - assemble global stiffness matrix
    assemble_global_stiffness();

    // step 3 Identify free DOFs (not FixedPin nodes)
    std::vector<int> free_dof_indices;
    free_dof_indices.reserve(total_dof);

    for (int i = 0; i < num_nodes; ++i)
    {
        // node dof indices
        int dof_x = i * 3;
        int dof_y = i * 3 + 1;
        int dof_theta = i * 3 + 2;

        if (nodes[i].constraint_type == Free)
        {
            // All 3 DOFs are free and must be solved for.
            free_dof_indices.push_back(dof_x);
            free_dof_indices.push_back(dof_y);
            free_dof_indices.push_back(dof_theta);
        }
        else if (nodes[i].constraint_type == Slider)
        {
            // x and y DOFs are free, theta is fixed
            free_dof_indices.push_back(dof_x);
            free_dof_indices.push_back(dof_y);
            free_dof_indices.push_back(dof_theta); // the beam can still rotate about the slider
        }
        else if (nodes[i].constraint_type == FixedPin)
        {
            // x and y DOFs are fixed, theta is free
            free_dof_indices.push_back(dof_theta);
            // dof_x and dof_y are constrained
        }
        else if (nodes[i].constraint_type == Fixed)
        {
            // All DOFs are fixed; do nothing
            continue;
        }
    }

    int num_free_dofs = static_cast<int>(free_dof_indices.size());

    if (num_free_dofs == 0)
    {
        std::cout << "No free DOFs to solve!" << std::endl;
        return -1;
    }

    // step 4 Create reduced stiffness matrix and force vector
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

    // step 5 Identify slider nodes
    std::vector<int> slider_nodes;
    for (int i = 0; i < num_nodes; ++i)
    {
        if (nodes[i].constraint_type == Slider)
            slider_nodes.push_back(i);
    }

    int num_constraints = static_cast<int>(slider_nodes.size());

    // step 5.1 if there are no constraints, solve directly
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
    // step 5.2 else build constraint matrix
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

        // step 6 Create augmented saddle point system

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

        // Debugging to determine the condition number if this is

        // stupid high we know the system is ill-conditioned

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(saddle_matrix);

        double cond = svd.singularValues()(0) / svd.singularValues().tail(1)(0);

        std::cout << "Saddle approx cond num: " << cond << std::endl;

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
        // The three DOFs for node i are at indices i*3, i*3+1, and i*3+2
        double u = displacement(i * 3);             // x-displacement
        double v = displacement(i * 3 + 1);         // y-displacement
        double theta_rad = displacement(i * 3 + 2); // Rotation (in radians) <-- NEW

        // Convert rotation to degrees for readability
        double theta_deg = theta_rad * 180.0 / M_PI;

        double total_disp_mag = std::sqrt(u * u + v * v);

        // Updated output includes rotation
        std::cout << "  Node " << i << ": u=" << u << " m, v=" << v << " m, theta=" << theta_deg << " deg (total disp=" << total_disp_mag << " m)\n";

        if (nodes[i].constraint_type == Slider)
        {
            // The stored angle in the node is the direction *along* the slider.
            double slider_angle_rad = nodes[i].constraint_angle * static_cast<double>(M_PI) / 180.0f;

            // Direction cosines of the slider track
            double c_track = std::cos(slider_angle_rad);
            double s_track = std::sin(slider_angle_rad);

            // 1. Movement along the slider track (dot product of displacement vector and track vector)
            double along_slider = u * c_track + v * s_track;

            // 2. Movement perpendicular to the slider track (dot product of displacement and normal vector)
            // Note: The normal vector cosines are ( -s_track, c_track )
            double perp_slider = -u * s_track + v * c_track;

            // The total displacement of the node should be explained as movement along the track.
            std::cout << "    Movement along track (" << nodes[i].constraint_angle << "°): " << along_slider << " m\n";

            // The perpendicular movement should be very close to zero if the constraint worked.
            std::cout << "    Movement perpendicular to track: " << perp_slider << " m (should be ~0)\n";
        }
    }

    // Calculate reaction forces and moments at supports
    // R = K*u - F. This 3N vector contains forces (Fx, Fy) and moments (Mz)
    Eigen::VectorXd reactions = global_k_matrix * displacement - forces;
    std::cout << "\nReaction Forces & Moments (N, Nm):\n";

    double total_reaction_x = 0.0;
    double total_reaction_y = 0.0;
    double total_reaction_m = 0.0; // NEW: Total Reaction Moment

    for (int i = 0; i < num_nodes; ++i)
    {
        // Only constrained nodes (Fixed, FixedPin, Slider) have reactions
        if (nodes[i].constraint_type != Free)
        {
            std::string constraint_str;
            if (nodes[i].constraint_type == Fixed)
                constraint_str = "Fixed";
            else if (nodes[i].constraint_type == FixedPin)
                constraint_str = "FixedPin";
            else
                constraint_str = "Slider"; // Handles the Slider constraint

            double rx = reactions(i * 3);     // x-Reaction Force
            double ry = reactions(i * 3 + 1); // y-Reaction Force
            double rm = reactions(i * 3 + 2); // NEW: z-Reaction Moment

            std::cout << "  Node " << i << " (" << constraint_str << "): Fx=" << rx << " N, Fy=" << ry << " N, Mz=" << rm << " Nm\n";

            total_reaction_x += rx;
            total_reaction_y += ry;
            total_reaction_m += rm; // Summing total reaction moment
        }
    }

    // Equilibrium check
    double total_applied_x = 0.0;
    double total_applied_y = 0.0;
    double total_applied_m = 0.0;

    for (int i = 0; i < num_nodes; ++i)
    {
        total_applied_x += forces(i * 3);
        total_applied_y += forces(i * 3 + 1);
        total_applied_m += forces(i * 3 + 2);
    }

    std::cout << "\nEquilibrium Check:\n";
    std::cout << "  Applied Fx = " << total_applied_x << " N\n";
    std::cout << "  Applied Fy = " << total_applied_y << " N\n";
    std::cout << "  Applied Mz = " << total_applied_m << " Nm\n";
    std::cout << "  Reaction Fx = " << total_reaction_x << " N\n";
    std::cout << "  Reaction Fy = " << total_reaction_y << " N\n";
    std::cout << "  Reaction Mz = " << total_reaction_m << " Nm\n";
    std::cout << "  Balance (Fx): " << (total_applied_x + total_reaction_x) << " N (should be ~0)\n";
    std::cout << "  Balance (Fy): " << (total_applied_y + total_reaction_y) << " N (should be ~0)\n";
    std::cout << "  Balance (Mz): " << (total_applied_m + total_reaction_m) << " Nm (should be ~0)\n";

    // Calculate internal forces and stresses in beams
    std::cout << "\nBeam Internal Forces (N, tension positive):\n";

    max_stress = -1e10f;
    min_stress = 1e10f;

    std::cout << "\nBeam Internal Forces and Moments (N, Nm):\n";
    int index = 0;
    for (auto &beam : beams) // Renamed 'spring' to 'beam' for clarity
    {
        int n1 = beam.nodes[0];
        int n2 = beam.nodes[1];

        // --- 1. Get Element Displacement Vector (6 DOFs) ---
        // The vector must be 6x1 in the order: {u1, v1, theta1, u2, v2, theta2}
        Eigen::VectorXd element_disp(6);
        element_disp << displacement(n1 * 3), // u1
            displacement(n1 * 3 + 1),         // v1
            displacement(n1 * 3 + 2),         // theta1
            displacement(n2 * 3),             // u2
            displacement(n2 * 3 + 1),         // v2
            displacement(n2 * 3 + 2);         // theta2

        // --- 2. Calculate Global Element End Forces/Moments (6x1 vector) ---
        // F_global = K_global * u_global
        Eigen::VectorXd global_end_forces = beam.k_matrix * element_disp; // 6x6 * 6x1 -> 6x1

        // --- 3. Setup Transformation Matrix (T) ---
        const Node &node1 = nodes[n1];
        const Node &node2 = nodes[n2];
        double dx = node2.position[0] - node1.position[0];
        double dy = node2.position[1] - node1.position[1];
        double L = std::sqrt(dx * dx + dy * dy);
        double c = dx / L;
        double s = dy / L;

        // Build the full 6x6 transformation matrix T
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(6, 6);
        // Block 1 (Node 1)
        T(0, 0) = c;
        T(0, 1) = s;
        T(1, 0) = -s;
        T(1, 1) = c;
        T(2, 2) = 1.0;
        // Block 2 (Node 2)
        T(3, 3) = c;
        T(3, 4) = s;
        T(4, 3) = -s;
        T(4, 4) = c;
        T(5, 5) = 1.0;

        // --- 4. Transform Forces to Local Coordinates ---
        // F_local = T * F_global
        // This gives the internal forces/moments in the order {P1, V1, M1, P2, V2, M2}
        Eigen::VectorXd local_end_forces = T * global_end_forces;

        // --- 5. Extract and Store Results ---
        // Axial Force (P): Must be taken from the end of the element (P2) or checked for sign.
        // Axial Force is typically constant along the beam, and should equal -P1 (local_end_forces(0)) or P2 (local_end_forces(3)).
        // P = Axial Force (Tension is positive)
        double P = local_end_forces(3); // Use P2 (index 3)

        // Shear Force (V)
        // double V1 = local_end_forces(1); // Shear at Node 1

        // Bending Moments (M)
        double M1 = local_end_forces(2);                     // Moment at Node 1 (local theta1 DOF)
        double M2 = local_end_forces(5);                     // Moment at Node 2 (local theta2 DOF)
        double max_M = std::max(std::abs(M1), std::abs(M2)); // Max moment magnitude

        // Update the beam object to store the results
        beam.axial_force = P;
        beam.max_moment = max_M;

        // Output Internal Forces and Moments
        std::cout << "  Beam " << index << " (nodes " << n1 << "-" << n2 << "): P=" << P << " N, M1=" << M1 << " Nm, M2=" << M2 << " Nm\n";

        // --- 6. Calculate Combined Stress ---
        const BeamProfile &shape = beam_profiles_list[beam.shape_idx];

        // Ensure Z_modulus is non-zero to avoid division by zero
        if (std::abs(shape.section_modulus) < 1e-12)
        {
            // Handle pure axial stress case or pure shear case if Z is effectively zero (e.g., truss I=0)
            beam.stress = P / shape.area;
        }
        else
        {
            // Combined Stress = Axial Stress (P/A) +/- Bending Stress (M/Z)
            // We calculate the maximum absolute stress using the max moment magnitude
            double axial_stress = P / shape.area;
            double bending_stress_max = max_M / shape.section_modulus;

            // Max tension stress (P/A + M/Z)
            double stress_tension = axial_stress + bending_stress_max;
            // Max compression stress (P/A - M/Z)
            double stress_compression = axial_stress - bending_stress_max;

            // Store the stress with correct sign: positive for tension, negative for compression
            if (std::abs(stress_tension) > std::abs(stress_compression))
                beam.stress = stress_tension;
            else
                beam.stress = stress_compression;
        }

        // Update global min/max stress trackers
        max_stress = std::max(max_stress, beam.stress);
        min_stress = std::min(min_stress, beam.stress);

        ++index;
    }

    std::cout << "\nBeam Maximum Absolute Combined Stresses (MPa):\n";
    index = 0;
    for (const auto &beam : beams)
    {
        std::cout << "  Beam " << index << " (nodes " << beam.nodes[0] << "-" << beam.nodes[1]
                  << "): " << beam.stress << " MPa\n";
        ++index;
    }

    std::cout << "\nStress Range: " << min_stress << " to " << max_stress << " MPa (Max Absolute Combined Stress)\n";

    return 0;
}

void FEMSystem::assemble_global_stiffness()
{
    int num_nodes = static_cast<int>(nodes.size());
    int total_dof = num_nodes * 3;
    global_k_matrix = Eigen::MatrixXd::Zero(total_dof, total_dof);

    // Assemble contributions from each spring
    for (const auto &spring : beams)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        int dof1_x = n1 * 3;
        int dof1_y = n1 * 3 + 1;
        int dof1_theta = n1 * 3 + 2;
        int dof2_x = n2 * 3;
        int dof2_y = n2 * 3 + 1;
        int dof2_theta = n2 * 3 + 2;

        int dofs[6] = {dof1_x, dof1_y, dof1_theta, dof2_x, dof2_y, dof2_theta};

        // Add element stiffness matrix to global matrix
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                global_k_matrix(dofs[i], dofs[j]) += spring.k_matrix(i, j);
            }
        }
    }

    std::cout << "Global Stiffness Matrix (" << total_dof << "x" << total_dof << "):\n"
              << global_k_matrix << std::endl;
}
