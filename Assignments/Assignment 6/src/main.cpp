#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/Graphics/Vertex.hpp>
#include <imgui.h>
#include <imgui-SFML.h>
#include <Eigen/Eigen>
#include <vector>
#include <algorithm>
#include <cstdint>

float scale_x = 1000;
float scale_y = 400;
float scales[2] = {scale_x, scale_y};

enum ConstraintType
{
    Free,
    Fixed,
    Slider
};

class Node
{
public:
    int id;
    float position[2]; // position in 2D space (x, y)
    ConstraintType constraint_type;
    float constraint_angle; // Angle in degrees (0 for horizontal, 90 for vertical)

    Node(int i, float x, float y, ConstraintType ct = Free, float angle = 0.0f) : id(i), constraint_type(ct), constraint_angle(angle)
    {
        position[0] = x;
        position[1] = y;
    }

    Eigen::Matrix2f get_transformation_matrix() const
    {
        Eigen::Matrix2f T = Eigen::Matrix2f::Identity();
        if (constraint_type == Slider)
        {
            float rad = constraint_angle * M_PI / 180.0f;
            float c = std::cos(rad);
            float s = std::sin(rad);
            T << c, s,
                -s, c;
        }
        return T;
    }
};

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
    void compute_stiffness(const std::vector<Node> &node_list)
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
};

class SpringSystem
{
public:
    std::vector<Node> nodes;
    std::vector<Spring> springs;
    Eigen::MatrixXd global_k_matrix;
    Eigen::VectorXd forces;       // forces in x and y directions [F1x, F1y, F2x, F2y, ...]
    Eigen::VectorXd displacement; // displacements in x and y [u1, v1, u2, v2, ...]
    float max_stress;
    float min_stress;

    SpringSystem(const std::vector<Node> &n, const std::vector<Spring> &s) : nodes(n), springs(s)
    {
        max_stress = 0.0f;
        min_stress = 0.0f;

        // Resize displacement and forces vectors
        int total_dof = nodes.size() * 2; // 2 DOF per node (x and y)
        displacement = Eigen::VectorXd::Zero(total_dof);
        forces = Eigen::VectorXd::Zero(total_dof);

        // step 1 - compute stiffness matrices for each spring
        // Compute spring stiffness matrices
        for (auto &spring : springs)
        {
            spring.compute_stiffness(nodes);
        }

        // step 2 - assemble global stiffness matrix
        assemble_global_stiffness();
    }

    // step 3 - attempt to solve the system
    int solve_system()
    {
        int num_nodes = nodes.size();
        int total_dof = num_nodes * 2; // 2 DOF per node (x and y)

        // Identify free DOFs (not fixed nodes)
        std::vector<int> free_dof_indices;

        for (int i = 0; i < num_nodes; ++i)
        {
            if (nodes[i].constraint_type != Fixed)
            {
                free_dof_indices.push_back(i * 2);     // x DOF
                free_dof_indices.push_back(i * 2 + 1); // y DOF
            }
        }

        int num_free_dofs = free_dof_indices.size();

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

        // Identify slider nodes and create constraint matrix
        std::vector<int> slider_nodes;
        for (int i = 0; i < num_nodes; ++i)
        {
            if (nodes[i].constraint_type == Slider)
                slider_nodes.push_back(i);
        }

        int num_constraints = slider_nodes.size();

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
            // Build constraint matrix in reduced coordinates
            Eigen::MatrixXd C_r = Eigen::MatrixXd::Zero(num_constraints, num_free_dofs);

            for (int ci = 0; ci < num_constraints; ++ci)
            {
                int node_id = slider_nodes[ci];
                const Node &node = nodes[node_id];

                float theta = node.constraint_angle * M_PI / 180.0f;
                float a_x = -std::sin(theta);
                float a_y = std::cos(theta);

                // Find the position of this node's DOFs in the reduced system
                for (int j = 0; j < num_free_dofs; ++j)
                {
                    int global_dof = free_dof_indices[j];
                    if (global_dof == node_id * 2)
                    {
                        C_r(ci, j) = a_x;
                    }
                    else if (global_dof == node_id * 2 + 1)
                    {
                        C_r(ci, j) = a_y;
                    }
                }
            }

            // Create augmented saddle point system
            // [K_r   C_r^T] [u]   [F_r]
            // [C_r    0   ] [λ] = [0  ]

            int augmented_size = num_free_dofs + num_constraints;
            Eigen::MatrixXd saddle_matrix = Eigen::MatrixXd::Zero(augmented_size, augmented_size);
            Eigen::VectorXd saddle_rhs = Eigen::VectorXd::Zero(augmented_size);

            // Fill in K_r block
            saddle_matrix.block(0, 0, num_free_dofs, num_free_dofs) = K_r;

            // Fill in C_r^T block (upper right)
            saddle_matrix.block(0, num_free_dofs, num_free_dofs, num_constraints) = C_r.transpose();

            // Fill in C_r block (lower left)
            saddle_matrix.block(num_free_dofs, 0, num_constraints, num_free_dofs) = C_r;

            // Fill in RHS
            saddle_rhs.head(num_free_dofs) = F_r;
            // Constraint RHS is zero (already initialized)

            // Solve the augmented system
            Eigen::VectorXd full_solution = saddle_matrix.fullPivLu().solve(saddle_rhs);

            // Extract displacements (first num_free_dofs entries)
            Eigen::VectorXd u_r = full_solution.head(num_free_dofs);

            // Extract Lagrange multipliers (last num_constraints entries)
            Eigen::VectorXd lagrange_multipliers = full_solution.tail(num_constraints);

            // Map back to global displacement vector
            displacement = Eigen::VectorXd::Zero(total_dof);
            for (int i = 0; i < num_free_dofs; ++i)
            {
                displacement(free_dof_indices[i]) = u_r(i);
            }

            // Print constraint violations for debugging
            Eigen::VectorXd constraint_check = C_r * u_r;
            std::cout << "\nConstraint violations (should be ~0):" << std::endl;
            for (int i = 0; i < num_constraints; ++i)
            {
                std::cout << "  Slider node " << slider_nodes[i] << ": " << constraint_check(i) << std::endl;
            }

            std::cout << "\nLagrange multipliers (constraint forces):" << std::endl;
            for (int i = 0; i < num_constraints; ++i)
            {
                std::cout << "  Slider node " << slider_nodes[i] << ": " << lagrange_multipliers(i) << " N" << std::endl;
            }
        }

        // Print solution
        std::cout << "\n=== SOLUTION ===" << std::endl;
        std::cout << "Displacements (m):" << std::endl;
        for (int i = 0; i < num_nodes; ++i)
        {
            std::cout << "  Node " << i << ": u=" << displacement(i * 2)
                      << " m, v=" << displacement(i * 2 + 1) << " m" << std::endl;
        }

        // Calculate reaction forces at fixed and slider nodes
        Eigen::VectorXd reactions = global_k_matrix * displacement - forces;
        std::cout << "\nReaction Forces (N):" << std::endl;
        for (int i = 0; i < num_nodes; ++i)
        {
            if (nodes[i].constraint_type != Free)
            {
                std::string constraint_str = (nodes[i].constraint_type == Fixed) ? "Fixed" : "Slider";
                std::cout << "  Node " << i << " (" << constraint_str << "): Fx=" << reactions(i * 2)
                          << " N, Fy=" << reactions(i * 2 + 1) << " N" << std::endl;
            }
        }

        // Calculate internal forces and stresses in springs
        std::cout << "\nSpring Internal Forces (N, tension positive):" << std::endl;
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

            float dx = nodes[n2].position[0] - nodes[n1].position[0];
            float dy = nodes[n2].position[1] - nodes[n1].position[1];
            float length = std::sqrt(dx * dx + dy * dy);
            float c = dx / length;
            float s = dy / length;

            // Axial force (negative of force at first node)
            float axial_force = -(c * element_forces(0) + s * element_forces(1));

            // Stress = Force / Area
            spring.stress = axial_force / spring.A / 1e6; // in MPa

            max_stress = std::max(max_stress, spring.stress);
            min_stress = std::min(min_stress, spring.stress);

            std::cout << "  Spring " << spring.id << " (nodes " << n1 << "-" << n2
                      << "): " << axial_force << " N" << std::endl;
        }

        // Print axial stresses in springs
        std::cout << "\nSpring Axial Stresses (MPa):" << std::endl;
        for (const auto &spring : springs)
        {
            std::cout << "  Spring " << spring.id << " (nodes " << spring.nodes[0] << "-" << spring.nodes[1]
                      << "): " << spring.stress << " MPa" << std::endl;
        }

        std::cout << "\nStress Range: " << min_stress << " to " << max_stress << " MPa" << std::endl;

        return 0;
    }

    void assemble_global_stiffness()
    {
        int num_nodes = nodes.size();
        int total_dof = num_nodes * 2;
        global_k_matrix = Eigen::MatrixXd::Zero(total_dof, total_dof);

        // Assemble contributions from each spring
        for (const auto &spring : springs)
        {
            int n1 = spring.nodes[0];
            int n2 = spring.nodes[1];

            // DOF indices: node i has DOFs [2*i, 2*i+1] for [x, y]
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

    // MPC (Multi-Point Constraint) for slider nodes
    void generate_constraint_row(Eigen::MatrixXd &C, int row_index, int node_id)
    {
        // get the node and its angle and we can compute the MCP
        const Node &node = nodes[node_id];
        if (node.constraint_type != Slider)
            return; // MPC only for slider nodes

        float theta = node.constraint_angle * M_PI / 180.0f;
        float a_x = -std::sin(theta);
        float a_y = std::cos(theta);

        int dof_x = node.id * 2;
        int dof_y = node.id * 2 + 1;

        // create the constraint row
        C(row_index, dof_x) = a_x;
        C(row_index, dof_y) = a_y;
    }
};

sf::Color getStressColor(float stress, float min_stress, float max_stress)
{
    float normalized;
    float abs_max = std::max(std::abs(min_stress), std::abs(max_stress));

    if (abs_max < 1e-6f)
    {
        return sf::Color::White;
    }

    normalized = stress / abs_max;

    if (normalized < 0)
    {
        float t = -normalized;
        uint8_t r = static_cast<uint8_t>(255 * (1 - t));
        uint8_t g = static_cast<uint8_t>(255 * (1 - t));
        uint8_t b = 255;
        return sf::Color(r, g, b);
    }
    else
    {
        float t = normalized;
        uint8_t r = 255;
        uint8_t g = static_cast<uint8_t>(255 * (1 - t));
        uint8_t b = static_cast<uint8_t>(255 * (1 - t));
        return sf::Color(r, g, b);
    }
}

void draw_system(sf::RenderWindow &window, const SpringSystem &system, float offset_x, float offset_y)
{
    // Draw springs with stress-based colors
    for (const auto &spring : system.springs)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        sf::Color springColor = getStressColor(spring.stress, system.min_stress, system.max_stress);

        sf::Vertex line[2];
        line[0].position = sf::Vector2f(offset_x + system.nodes[n1].position[0] * scale_x + system.displacement(n1 * 2) * scale_x,
                                        offset_y - system.nodes[n1].position[1] * scale_y - system.displacement(n1 * 2 + 1) * scale_y);
        line[0].color = springColor;

        line[1].position = sf::Vector2f(offset_x + system.nodes[n2].position[0] * scale_x + system.displacement(n2 * 2) * scale_x,
                                        offset_y - system.nodes[n2].position[1] * scale_y - system.displacement(n2 * 2 + 1) * scale_y);
        line[1].color = springColor;

        window.draw(line, 2, sf::PrimitiveType::Lines);
    }

    // Draw nodes
    for (const auto &node : system.nodes)
    {
        if (node.constraint_type == Fixed)
        {
            sf::RectangleShape square(sf::Vector2f(10, 10));
            square.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x - 5,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y - 5));
            square.setFillColor(sf::Color::Red);
            window.draw(square);
        }

        else if (node.constraint_type == Slider)
        {
            sf::ConvexShape triangle(3); // Create a triangle with 3 points
            triangle.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y));

            // Direction perpendicular to the slider (add 90° to the constraint angle)
            float perp_angle_rad = (node.constraint_angle + 90.0f) * M_PI / 180.0f;
            float dx = std::cos(perp_angle_rad);
            float dy = -std::sin(perp_angle_rad); // Negative because y-axis is inverted in SFML

            // Triangle pointing in the perpendicular direction
            triangle.setPoint(0, sf::Vector2f(0, 0));                                 // Tip of the triangle
            triangle.setPoint(1, sf::Vector2f(-dx * 10 - dy * 5, -dy * 10 + dx * 5)); // Base left
            triangle.setPoint(2, sf::Vector2f(-dx * 10 + dy * 5, -dy * 10 - dx * 5)); // Base right
            triangle.setFillColor(sf::Color::Yellow);
            window.draw(triangle);
        }
        else
        {
            sf::CircleShape circle(5);
            circle.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x - 5,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y - 5));
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);
        }
    }

    // Draw forces as arrows on nodes
    for (int i = 0; i < system.nodes.size(); ++i)
    {
        const Node &node = system.nodes[i];
        float fx = system.forces(i * 2) / 1000.0f;      // Scale down for visualization
        float fy = -system.forces(i * 2 + 1) / 1000.0f; // Invert y for SFML
        if (std::abs(fx) < 1e-3f && std::abs(fy) < 1e-3f)
            continue; // Skip zero forces
        sf::Vector2f start(offset_x + node.position[0] * scale_x + system.displacement(i * 2) * scale_x,
                           offset_y - node.position[1] * scale_y - system.displacement(i * 2 + 1) * scale_y);
        sf::Vector2f end = start + sf::Vector2f(fx, fy);
        sf::Vertex line[2];
        line[0].position = start;
        line[0].color = sf::Color::Magenta;
        line[1].position = end;
        line[1].color = sf::Color::Magenta;

        // Draw arrowhead
        sf::Vector2f direction = end - start;
        float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
        if (length > 0)
        {
            sf::Vector2f unit_dir = direction / length;
            sf::Vector2f perp_dir(-unit_dir.y, unit_dir.x);
            float arrow_size = 5.0f;
            sf::Vertex arrow[3];
            arrow[0].position = end;
            arrow[0].color = sf::Color::Magenta;
            arrow[1].position = end - unit_dir * arrow_size + perp_dir * (arrow_size / 2);
            arrow[1].color = sf::Color::Magenta;
            arrow[2].position = end - unit_dir * arrow_size - perp_dir * (arrow_size / 2);
            arrow[2].color = sf::Color::Magenta;
            window.draw(arrow, 3, sf::PrimitiveType::Triangles);
        }

        window.draw(line, 2, sf::PrimitiveType::Lines);
    }
}

int main()
{
    sf::Vector2u window_size(1200, 600);

    sf::RenderWindow window(sf::VideoMode(window_size), "2D Spring System Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    sf::Clock deltaClock;

    // Define nodes
    std::vector<Node> nodes = {
        Node(0, 0.0f, 0.0f, Fixed),                  // Fixed node at origin
        Node(1, 0.0f, 1.0f, Slider, 360.0f - 90.0f), // Slider node
        Node(2, 0.5f, 1.0f, Slider, 60.0f),          // Slider node
    };

    const double A = 6.0e-4; // Cross-sectional area in m^2
    const double E = 210e9;  // Young's modulus in Pa

    // Define springs
    std::vector<Spring> springs = {
        Spring(0, 0, 1, A, E),
        Spring(1, 1, 2, A, E),
        Spring(2, 2, 0, A, E),
    };

    SpringSystem spring_system(nodes, springs);

    // Apply a 400 kN downward force at node 2
    spring_system.forces(2 * 2 + 1) = -400000.0f; // Fy at node 2 (negative = downward)

    std::cout << "\nApplied Forces (N):" << std::endl;
    for (int i = 0; i < nodes.size(); ++i)
    {
        std::cout << "  Node " << i << ": Fx=" << spring_system.forces(i * 2)
                  << " N, Fy=" << spring_system.forces(i * 2 + 1) << " N" << std::endl;
    }

    spring_system.solve_system();

    while (window.isOpen())
    {
        while (const std::optional<sf::Event> event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }

            if (const auto *keyPressed = event->getIf<sf::Event::KeyPressed>())
            {
                if (keyPressed->code == sf::Keyboard::Key::Escape)
                {
                    window.close();
                }
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        ImGui::Begin("System Controls");

        ImGui::Text("2D Truss System Solver");
        ImGui::Separator();

        // Scaling controls
        ImGui::SliderFloat("Scale X", &scale_x, 50.0f, 2000.0f);
        ImGui::SliderFloat("Scale Y", &scale_y, 50.0f, 2000.0f);

        ImGui::Separator();

        // Force controls for non-fixed nodes
        ImGui::Text("Apply Forces (N):");
        bool forces_changed = false;
        for (int i = 0; i < spring_system.nodes.size(); ++i)
        {
            if ((spring_system.nodes[i].constraint_type) == Free)
            {
                float fx = spring_system.forces(i * 2);
                float fy = spring_system.forces(i * 2 + 1);
                ImGui::PushID(i);

                // Fx and Fy controls
                ImGui::Text("Node %d Forces:", i);

                if (ImGui::SliderFloat("Fx", &fx, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                if (ImGui::SliderFloat("Fy", &fy, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }

                ImGui::Separator();
                ImGui::PopID();
            }
            else if (spring_system.nodes[i].constraint_type == Slider)
            {
                float fx = spring_system.forces(i * 2);
                float fy = spring_system.forces(i * 2 + 1);

                // Forces are already in global coordinates, no transformation needed

                ImGui::PushID(i);

                // Fx and Fy controls
                ImGui::Text("Node %d Forces (Slider):", i);

                ImGui::Columns(2, nullptr, false); // Create two columns for sliders and inputs

                ImGui::SetColumnWidth(0, ImGui::GetWindowWidth() * 0.8f); // Sliders take up most of the space
                ImGui::Text("Fx Slider:");
                if (ImGui::SliderFloat("##FxSlider", &fx, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                ImGui::NextColumn();
                ImGui::Text("Fx Input:");
                if (ImGui::InputFloat("##FxInput", &fx, 0.0f, 0.0f, "%.1f"))
                {
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                ImGui::NextColumn();

                ImGui::SetColumnWidth(0, ImGui::GetWindowWidth() * 0.8f); // Sliders take up most of the space
                ImGui::Text("Fy Slider:");
                if (ImGui::SliderFloat("##FySlider", &fy, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }
                ImGui::NextColumn();
                ImGui::Text("Fy Input:");
                if (ImGui::InputFloat("##FyInput", &fy, 0.0f, 0.0f, "%.1f"))
                {
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }
                ImGui::Columns(1); // Reset to single column layout

                ImGui::Separator();
                ImGui::PopID();
            }
        }

        if (forces_changed)
        {
            spring_system.solve_system();
        }

        ImGui::Separator();

        // Display solution info
        ImGui::Text("Solution:");
        for (int i = 0; i < spring_system.nodes.size(); ++i)
        {
            ImGui::Text("Node %d: u=%.6f m, v=%.6f m", i,
                        spring_system.displacement(i * 2),
                        spring_system.displacement(i * 2 + 1));
        }

        ImGui::Separator();

        // Display spring stresses with color legend
        ImGui::Text("Spring Stresses (MPa):");
        ImGui::Text("Blue = Compression, Red = Tension");
        for (const auto &spring : spring_system.springs)
        {
            sf::Color color = getStressColor(spring.stress, spring_system.min_stress, spring_system.max_stress);
            ImGui::TextColored(ImVec4(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f, 1.0f),
                               "Spring %d: %.2f MPa", spring.id, spring.stress);
        }

        ImGui::Separator();
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                    1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        ImGui::End();

        window.clear(sf::Color(30, 30, 30));

        static float origin_x = 100.0f;
        static float origin_y = 300.0f;

        // Draw the system
        draw_system(window, spring_system, origin_x, origin_y);

        // let the user drag the origin with the mouse, only if clicking near the origin
        if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
        {
            sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
            float distance_to_origin = std::sqrt(std::pow(mouse_pos.x - origin_x, 2) + std::pow(mouse_pos.y - origin_y, 2));

            // Allow dragging only if the mouse is within a certain radius of the origin
            const float drag_radius = 20.0f; // Radius in pixels
            if (distance_to_origin <= drag_radius)
            {
                origin_x = static_cast<float>(mouse_pos.x);
                origin_y = static_cast<float>(mouse_pos.y);
            }
        }

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
