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

class Node
{
public:
    int id;
    float position[2]; // position in 2D space (x, y)
    bool is_fixed;

    Node(int i, float x, float y, bool fixed) : id(i), is_fixed(fixed)
    {
        position[0] = x;
        position[1] = y;
    }
};

class Spring
{
public:
    Spring(int i, int n1, int n2) : id(i)
    {
        nodes[0] = n1;
        nodes[1] = n2;
        stress = 0.0f;
    }

    int id;
    int nodes[2]; // node indices
    double k;     // stiffness (EA/L)
    float stress; // axial stress in MPa

    Eigen::Matrix4f k_matrix; // 4x4 stiffness matrix in global coordinates
                              // [u1, v1, u2, v2]

    void compute_stiffness(const std::vector<Node> &node_list)
    {
        int n1 = nodes[0];
        int n2 = nodes[1];

        float dx = node_list[n2].position[0] - node_list[n1].position[0];
        float dy = node_list[n2].position[1] - node_list[n1].position[1];
        float length = std::sqrt(dx * dx + dy * dy);

        // EA/L where E=210 GPa, A=6e-4 m^2
        k = (6.0e-4 * 210e9) / length;

        // Direction cosines
        float c = dx / length; // cos(theta)
        float s = dy / length; // sin(theta)

        // local stiffness matrix for a 2D truss element
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
    Eigen::MatrixXf global_k_matrix;
    Eigen::VectorXf forces;       // forces in x and y directions [F1x, F1y, F2x, F2y, ...]
    Eigen::VectorXf displacement; // displacements in x and y [u1, v1, u2, v2, ...]
    float max_stress;
    float min_stress;

    SpringSystem(const std::vector<Node> &n, const std::vector<Spring> &s) : nodes(n), springs(s)
    {
        max_stress = 0.0f;
        min_stress = 0.0f;

        // Compute spring stiffness matrices
        for (auto &spring : springs)
        {
            spring.compute_stiffness(nodes);
        }

        assemble_global_stiffness();
    }

    int solve_system()
    {
        int num_nodes = nodes.size();
        int total_dof = num_nodes * 2; // 2 DOF per node (x and y)

        // Identify free DOFs
        std::vector<int> free_dof_indices;
        for (int i = 0; i < num_nodes; ++i)
        {
            if (!nodes[i].is_fixed)
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

        // Create reduced system
        Eigen::MatrixXf reduced_k_matrix = Eigen::MatrixXf::Zero(num_free_dofs, num_free_dofs);
        Eigen::VectorXf reduced_forces = Eigen::VectorXf::Zero(num_free_dofs);

        for (int i = 0; i < num_free_dofs; ++i)
        {
            for (int j = 0; j < num_free_dofs; ++j)
                reduced_k_matrix(i, j) = global_k_matrix(free_dof_indices[i], free_dof_indices[j]);
            reduced_forces(i) = forces(free_dof_indices[i]);
        }

        std::cout << "\nFree DOF indices: ";
        for (int idx : free_dof_indices)
            std::cout << idx << " ";
        std::cout << "\n\nReduced K matrix:\n"
                  << reduced_k_matrix << std::endl;
        std::cout << "\nReduced forces:\n"
                  << reduced_forces << std::endl;

        // Solve K * u = F
        Eigen::VectorXf reduced_displacements = reduced_k_matrix.colPivHouseholderQr().solve(reduced_forces);

        // Map back to global displacement vector
        displacement = Eigen::VectorXf::Zero(total_dof);
        for (int i = 0; i < num_free_dofs; ++i)
            displacement(free_dof_indices[i]) = reduced_displacements(i);

        std::cout << "\n=== SOLUTION ===" << std::endl;
        std::cout << "Displacements (m):" << std::endl;
        for (int i = 0; i < num_nodes; ++i)
        {
            std::cout << "  Node " << i << ": u=" << displacement(i * 2)
                      << " m, v=" << displacement(i * 2 + 1) << " m" << std::endl;
        }

        // Calculate reaction forces at fixed nodes
        Eigen::VectorXf reactions = global_k_matrix * displacement - forces;
        std::cout << "\nReaction Forces (N):" << std::endl;
        for (int i = 0; i < num_nodes; ++i)
        {
            if (nodes[i].is_fixed)
            {
                std::cout << "  Node " << i << ": Fx=" << reactions(i * 2)
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

            Eigen::Vector4f element_disp;
            element_disp << displacement(n1 * 2), displacement(n1 * 2 + 1),
                displacement(n2 * 2), displacement(n2 * 2 + 1);

            Eigen::Vector4f element_forces = spring.k_matrix * element_disp;

            float dx = nodes[n2].position[0] - nodes[n1].position[0];
            float dy = nodes[n2].position[1] - nodes[n1].position[1];
            float length = std::sqrt(dx * dx + dy * dy);
            float c = dx / length;
            float s = dy / length;

            // Axial force (negative of force at first node)
            float axial_force = -(c * element_forces(0) + s * element_forces(1));

            // Stress = Force / Area, Area = 6e-4 m^2
            spring.stress = axial_force / 6.0e-4 / 1e6; // in MPa

            max_stress = std::max(max_stress, spring.stress);
            min_stress = std::min(min_stress, spring.stress);

            std::cout << "  Spring " << spring.id << " (nodes " << n1 << "-" << n2
                      << "): " << axial_force << " N" << std::endl;
        }

        // calculate axial stresses in springs
        std::cout << "\nSpring Axial Stresses (MPa):" << std::endl;
        for (const auto &spring : springs)
        {
            std::cout << "  Spring " << spring.id << " (nodes " << spring.nodes[0] << "-" << spring.nodes[1]
                      << "): " << spring.stress << " MPa" << std::endl;
        }

        std::cout << "\nStress Range: " << min_stress << " to " << max_stress << " MPa" << std::endl;

        return 0;
    }

private:
    void assemble_global_stiffness()
    {
        int num_nodes = nodes.size();
        int total_dof = num_nodes * 2;
        global_k_matrix = Eigen::MatrixXf::Zero(total_dof, total_dof);

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

        forces = Eigen::VectorXf::Zero(total_dof);
        displacement = Eigen::VectorXf::Zero(total_dof);

        std::cout << "Global Stiffness Matrix (" << total_dof << "x" << total_dof << "):\n"
                  << global_k_matrix << std::endl;
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
        sf::CircleShape circle(5);
        circle.setPosition(sf::Vector2f(
            offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x - 5,
            offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y - 5));
        circle.setFillColor(node.is_fixed ? sf::Color::Red : sf::Color::Blue);
        window.draw(circle);
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
        Node(0, 0.0f, 0.0f, true),      // Fixed node at origin
        Node(1, 1.0f, 0.0f, false),     // Free node
        Node(2, 1.0f, 0.57735f, false), // Free node
        Node(3, 0.0f, 0.57735f, true),  // Fixed node
    };

    // Define springs
    std::vector<Spring> springs = {
        Spring(0, 0, 1),
        Spring(1, 1, 2),
        Spring(2, 2, 3),
        Spring(3, 0, 2),
    };

    SpringSystem spring_system(nodes, springs);

    // Apply a force of 500kN downward at node 1 and 400kN rightward at node 2
    spring_system.forces(1 * 2 + 1) = -500000.0; // Fy at node 1
    spring_system.forces(2 * 2) = 400000.0;      // Fx at node 2

    std::cout << "\nApplied Forces (N):" << std::endl;
    for (int i = 0; i < nodes.size(); ++i)
    {
        std::cout << "  Node " << i << ": Fx=" << spring_system.forces(i * 2)
                  << " N, Fy=" << spring_system.forces(i * 2 + 1) << " N" << std::endl;
    }

    for (const auto &spring : springs)
    {
        std::cout << "Spring " << spring.id << ": nodes "
                  << spring.nodes[0] << " → " << spring.nodes[1] << std::endl;
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
            if (!spring_system.nodes[i].is_fixed)
            {
                float fx = spring_system.forces(i * 2);
                float fy = spring_system.forces(i * 2 + 1);
                ImGui::PushID(i);

                // Fx controls
                ImGui::Text("Node %d Fx:", i);
                if (ImGui::SliderFloat("##FxSlider", &fx, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::InputFloat("##FxInput", &fx, 1000.0f, 10000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Fx +"))
                {
                    fx += 1000.0f;
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Fx -"))
                {
                    fx -= 1000.0f;
                    spring_system.forces(i * 2) = fx;
                    forces_changed = true;
                }

                // Fy controls
                ImGui::Text("Node %d Fy:", i);
                if (ImGui::SliderFloat("##FySlider", &fy, -1000000.0f, 1000000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::InputFloat("##FyInput", &fy, 1000.0f, 10000.0f, "%.1f"))
                {
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Fy +"))
                {
                    fy += 1000.0f;
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Fy -"))
                {
                    fy -= 1000.0f;
                    spring_system.forces(i * 2 + 1) = fy;
                    forces_changed = true;
                }

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

        // Draw the system
        draw_system(window, spring_system, 100, 300);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}