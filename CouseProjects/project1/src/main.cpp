#include <iostream>
#include <cmath>
#include <Eigen/Eigen>
#include <vector>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <imgui-SFML.h>
#include <imgui.h>

// Constants for scaling
float scale_x = 1000;
float scale_y = 4;
float scales[2] = {scale_x, scale_y};

// Class representing a node in the spring system
class Node
{
public:
    int id;
    float position;
    bool is_fixed;

    Node(int i, float pos, bool fixed) : id(i), position(pos), is_fixed(fixed) {}
};

// Class representing a spring in the spring system
class Spring
{
public:
    Spring(int i, int n1, int n2, double stiffness) : id(i), k(stiffness)
    {
        nodes[0] = n1;
        nodes[1] = n2;
        compute_stiffness();
    }

    int id;
    int nodes[2];             // node indices
    double k;                 // stiffness
    Eigen::Matrix2f k_matrix; // stiffness matrix

private:
    void compute_stiffness() { k_matrix << k, -k, -k, k; }
};

// Function to compute springyness
double compute_springyness(double E, double A, double L) { return (E * A) / L; }

// Function to compute width
double compute_width(double x)
{
    return 4.0 - x * (1.0 / 5.0); // linear taper from 4 to 2 over 10 units
}

// Class representing the spring system
class SpringSystem
{
public:
    std::vector<Node> nodes;
    std::vector<Spring> springs;
    Eigen::MatrixXf global_k_matrix;
    Eigen::VectorXf forces;
    Eigen::VectorXf displacements;

    SpringSystem(const std::vector<Node> &n, const std::vector<Spring> &s) : nodes(n), springs(s) { assemble_global_stiffness(); }

    int solve_system()
    {
        // Check to see if we have any fixed nodes
        int num_nodes = nodes.size();
        Eigen::VectorXf reduced_forces;
        Eigen::MatrixXf reduced_k_matrix;
        std::vector<int> free_node_indices;

        for (int i = 0; i < num_nodes; ++i)
        {
            if (!nodes[i].is_fixed)
            {
                free_node_indices.push_back(i);
            }
        }

        int num_free_nodes = free_node_indices.size();
        reduced_k_matrix = Eigen::MatrixXf::Zero(num_free_nodes, num_free_nodes);
        reduced_forces = Eigen::VectorXf::Zero(num_free_nodes);

        for (int i = 0; i < num_free_nodes; ++i)
        {
            for (int j = 0; j < num_free_nodes; ++j)
            {
                reduced_k_matrix(i, j) = global_k_matrix(free_node_indices[i], free_node_indices[j]);
            }
            reduced_forces(i) = forces(free_node_indices[i]);
        }

        // Solve for displacements
        Eigen::VectorXf reduced_displacements = reduced_k_matrix.colPivHouseholderQr().solve(reduced_forces);

        for (int i = 0; i < num_free_nodes; ++i)
        {
            displacements(free_node_indices[i]) = reduced_displacements(i);
        }

        // Solve for forces
        forces = global_k_matrix * displacements;

        return 0;
    }

private:
    void assemble_global_stiffness()
    {
        int num_nodes = nodes.size();
        global_k_matrix = Eigen::MatrixXf::Zero(num_nodes, num_nodes);
        for (const auto &spring : springs)
        {
            int n1 = spring.nodes[0];
            int n2 = spring.nodes[1];
            global_k_matrix(n1, n1) += spring.k_matrix(0, 0);
            global_k_matrix(n1, n2) += spring.k_matrix(0, 1);
            global_k_matrix(n2, n1) += spring.k_matrix(1, 0);
            global_k_matrix(n2, n2) += spring.k_matrix(1, 1);
        }
        forces = Eigen::VectorXf::Zero(num_nodes);
        displacements = Eigen::VectorXf::Zero(num_nodes);
    }
};

// Function to compute diameter
double compute_diameter(double x, double D, double d, double L) { return D - (D - d) * (x / L); }

// Function to compute area
double compute_area(double diameter) { return M_PI * std::pow(diameter / 2.0, 2); }

// Function to compute stress
double compute_stress(double force, double area) { return force / area; }

int main()
{
    // Constants
    const double E = 260e9; // Young's modulus in Pascals (N/m^2)
    const double P = 500;   // Applied load in Newtons
    const double D = 0.75;  // Major diameter in meters
    const double d = 0.25;  // Minor diameter in meters
    const double L = 1.0;   // Length in meters
    const int segments = 300;

    std::vector<Node> nodes;
    std::vector<Spring> springs;

    // Create nodes and springs
    for (int i = 1; i < segments; ++i)
    {
        float position = (L / segments) * i;
        bool is_fixed = (i == 1);
        nodes.emplace_back(i, position, is_fixed);

        double x = position - (L / segments);
        double width = compute_width(x);
        double diameter = compute_diameter(x, D, d, L);
        double area = compute_area(diameter);
        double springyness = compute_springyness(E, area, L / segments);
        springs.emplace_back(i - 1, i - 1, i, springyness);
    }

    nodes.emplace_back(segments, L, true); // Fixed node at the end

    std::cout << "Total nodes created: " << nodes.size() + 1 << std::endl; // +1 for the last fixed node
    std::cout << "Total springs created: " << springs.size() << std::endl; // add the last node

    SpringSystem spring_system(nodes, springs);

    float applied_force_location = 0;  // we want to set this to where the diff of the stress is minimized
    float stress_diff_at_location = 0; // apply the force and sweep it across the free nodes one at a time to find the change in reaction forces

    std::vector<float> stress_diffs;

    // Progress bar variables
    const int progress_bar_width = 50;
    std::cout << "Sweeping force application across nodes:\n[";

    for (int i = 0; i < segments; ++i)
    {
        if (!nodes[i].is_fixed)
        {
            spring_system.forces(i) = P;
            spring_system.solve_system();

            int start_node = 0;
            int end_node = segments - 1;
            float reaction_forces[2];
            reaction_forces[0] = spring_system.forces(start_node);
            reaction_forces[1] = spring_system.forces(end_node);

            float stress_diff = std::abs(compute_stress(reaction_forces[0], compute_area(compute_diameter(0.0, D, d, L))) - compute_stress(reaction_forces[1], compute_area(compute_diameter(L, D, d, L))));
            stress_diffs.push_back(stress_diff);
            if (i == 1 || stress_diff < stress_diff_at_location)
            {
                applied_force_location = i;
                stress_diff_at_location = stress_diff;
            }

            spring_system.displacements.setZero(); // reset displacements
            spring_system.forces.setZero();        // reset forces
        }

        // Progress bar update
        if (i % (segments / progress_bar_width) == 0)
        {
            std::cout << "#";
            std::cout.flush();
        }
    }

    std::cout << "] Done.\n";

    // Now we print the distance where the force should be applied and the stress difference at that location
    std::cout << "Optimal force application location: Node " << applied_force_location << " (Position: " << nodes[static_cast<int>(applied_force_location)].position << " meters)" << std::endl;
    std::cout << "Minimum stress difference at this location: " << stress_diff_at_location << " Pa" << std::endl;

    // use sfml to plot the stress diffs
    // Create a window
    sf::Vector2u window_size(1200, 400);
    sf::RenderWindow window(sf::VideoMode(window_size), "Stress Difference Plot");
    window.setFramerateLimit(60);
    // Initialize ImGui-SFML
    ImGui::SFML::Init(window);
    // Main loop
    while (window.isOpen())
    {
        while (const std::optional<sf::Event> event = window.pollEvent())
        {
            ImGui::SFML::ProcessEvent(window, *event);
            // Handle close event
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }

            // Handle other application-specific events
            if (const auto *keyPressed = event->getIf<sf::Event::KeyPressed>())
            {
                if (keyPressed->code == sf::Keyboard::Key::Escape)
                {
                    window.close();
                }
            }
        }
        ImGui::SFML::Update(window, sf::seconds(1.f / 60.f));
        window.clear(sf::Color::White);
        // Create a window for the plot
        // Create a window for the plot
        ImGui::Begin("Stress Difference Plot", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

        // Prepare data for plotting: x = position (meters), y = stress difference (Pa)
        static std::vector<float> x_positions;
        if (x_positions.size() != stress_diffs.size())
        {
            x_positions.resize(stress_diffs.size());
            for (size_t i = 0; i < stress_diffs.size(); ++i)
            {
                x_positions[i] = nodes[i].position;
            }
        }

        // Find min/max for scaling
        float min_stress = *std::min_element(stress_diffs.begin(), stress_diffs.end());
        float max_stress = *std::max_element(stress_diffs.begin(), stress_diffs.end());
        float min_pos = x_positions.front();
        float max_pos = x_positions.back();

        // Plot area
        ImVec2 plot_size(800, 400);
        ImVec2 plot_pos = ImGui::GetCursorScreenPos();
        ImDrawList *draw_list = ImGui::GetWindowDrawList();

        // Draw background
        draw_list->AddRectFilled(plot_pos, ImVec2(plot_pos.x + plot_size.x, plot_pos.y + plot_size.y),
                                 IM_COL32(40, 40, 40, 255));

        // Draw grid lines
        for (int i = 0; i <= 10; ++i)
        {
            float y = plot_pos.y + (plot_size.y / 10.0f) * i;
            draw_list->AddLine(ImVec2(plot_pos.x, y), ImVec2(plot_pos.x + plot_size.x, y),
                               IM_COL32(80, 80, 80, 255));
        }

        // Plot the data
        for (size_t i = 1; i < stress_diffs.size(); ++i)
        {
            // Map position to X coordinate
            float x1 = plot_pos.x + ((x_positions[i - 1] - min_pos) / (max_pos - min_pos)) * plot_size.x;
            float x2 = plot_pos.x + ((x_positions[i] - min_pos) / (max_pos - min_pos)) * plot_size.x;

            // Map stress to Y coordinate (inverted because screen Y increases downward)
            float y1 = plot_pos.y + plot_size.y - ((stress_diffs[i - 1] - min_stress) / (max_stress - min_stress)) * plot_size.y;
            float y2 = plot_pos.y + plot_size.y - ((stress_diffs[i] - min_stress) / (max_stress - min_stress)) * plot_size.y;

            draw_list->AddLine(ImVec2(x1, y1), ImVec2(x2, y2), IM_COL32(0, 255, 0, 255), 2.0f);
        }

        // Create invisible button for hover detection
        ImGui::InvisibleButton("plot", plot_size);

        // Show hovered value
        if (ImGui::IsItemHovered())
        {
            ImGuiIO &io = ImGui::GetIO();
            float mouse_x = io.MousePos.x - plot_pos.x;
            float mouse_y = io.MousePos.y - plot_pos.y;

            // Map mouse X back to position
            float hovered_pos = min_pos + (mouse_x / plot_size.x) * (max_pos - min_pos);

            // Find closest data point
            int closest_idx = 0;
            float min_dist = std::abs(x_positions[0] - hovered_pos);
            for (size_t i = 1; i < x_positions.size(); ++i)
            {
                float dist = std::abs(x_positions[i] - hovered_pos);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    closest_idx = i;
                }
            }

            // Draw tooltip
            ImGui::BeginTooltip();
            ImGui::Text("Position: %.4f meters", x_positions[closest_idx]);
            ImGui::Text("Stress Difference: %.2f Pa", stress_diffs[closest_idx]);
            ImGui::EndTooltip();

            // Draw crosshair at hovered point
            float px = plot_pos.x + ((x_positions[closest_idx] - min_pos) / (max_pos - min_pos)) * plot_size.x;
            float py = plot_pos.y + plot_size.y - ((stress_diffs[closest_idx] - min_stress) / (max_stress - min_stress)) * plot_size.y;
            draw_list->AddCircleFilled(ImVec2(px, py), 5.0f, IM_COL32(255, 255, 0, 255));
        }

        // Draw axis labels
        ImGui::Text("Position: %.3f m to %.3f m", min_pos, max_pos);
        ImGui::Text("Stress Difference: %.2f Pa to %.2f Pa", min_stress, max_stress);
        ImGui::Text("Optimal location: %.4f m (Stress diff: %.2f Pa)",
                    nodes[applied_force_location].position, stress_diff_at_location);

        ImGui::End();
        // Render ImGui
        ImGui::SFML::Render(window);
        window.display();
    }
    ImGui::SFML::Shutdown();

    return 0;
}