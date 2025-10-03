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

float scale_x = 1000;
float scale_y = 4;
float scales[2] = {scale_x, scale_y};

class Node
{
public:
    int id;
    float position;
    bool is_fixed;

    Node(int i, float pos, bool fixed) : id(i), position(pos), is_fixed(fixed) {}

    void draw(sf::RenderWindow &window) const
    {
        sf::CircleShape circle(10);
        circle.setFillColor(is_fixed ? sf::Color::Red : sf::Color::Green);
        circle.setPosition(sf::Vector2f(position * (scale_x / 1000.0f), static_cast<float>(window.getSize().y) / 2.0f - 10.0f));
        window.draw(circle);
    }
};

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
    int nodes[2]; // node indices
    double k;     // stiffness

    Eigen::Matrix2f k_matrix; // stiffness matrix

private:
    void compute_stiffness()
    {
        k_matrix << k, -k,
            -k, k;
    }
};

double compute_springyness(double E, double A, double L)
{
    return (E * A) / L;
}

double compute_width(double x)
{
    return 4.0 - x * (1.0 / 5.0); // linear taper from 4 to 2 over 10 units
}

class SpringSystem
{
public:
    std::vector<Node> nodes;
    std::vector<Spring> springs;
    Eigen::MatrixXf global_k_matrix;
    Eigen::VectorXf forces;
    Eigen::VectorXf displacements;

    // the user defines the spring system by providing a vector of the nodes and a vector of the springs
    SpringSystem(const std::vector<Node> &n, const std::vector<Spring> &s) : nodes(n), springs(s)
    {
        assemble_global_stiffness();
    }

    int solve_system()
    {
        // check to see if we have any fixed nodes
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

        // solve for forces
        forces = global_k_matrix * displacements;

        // std::cout << "Displacements:\n"
        //           << displacements << std::endl;

        // std::cout << "Forces:\n"
        //           << forces << std::endl;
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

        // std::cout << "Global Stiffness Matrix:\n"
        //           << global_k_matrix << std::endl;
    }
};

double compute_diameter(double x, double D, double d, double L)
{
    return D - (D - d) * (x / L);
}

double compute_area(double diameter)
{
    return M_PI * std::pow(diameter / 2.0, 2);
}

double compute_stress(double force, double area)
{
    return force / area;
}

int main()
{
    sf::Vector2u window_size(1200, 400);

    const double E = 260e9; // Young's modulus in Pascals (N/m^2)
    const double P = 500;   // Applied load in Newtons
    const double D = 0.75;  // Major diameter in meters
    const double d = 0.25;  // Minor diameter in meters
    const double L = 1.0;   // Length in meters

    const int segments = 1000;

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

        // std::cout << "Node " << i << ": Position = " << position << ", Fixed = " << is_fixed << std::endl;
    }

    std::cout << "Total nodes created: " << nodes.size() + 1 << std::endl; // +1 for the last fixed node
    std::cout << "Total springs created: " << springs.size() << std::endl;

    // add the last node
    nodes.emplace_back(segments, L, true); // Fixed node at the end
    std::cout << "Node " << segments << ": Position = " << L << ", Fixed = " << 1 << std::endl;

    SpringSystem spring_system(nodes, springs);

    float applied_force_location = 0; // we want to set this to where the diff of the stress is minimized
    float stress_diff_at_location = 0;

    // apply the force and sweep it across the free nodes one at a time to find the change in reaction forces
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
            float stress_diff = std::abs(compute_stress(reaction_forces[0], compute_area(compute_diameter(0, D, d, L))) -
                                         compute_stress(reaction_forces[1], compute_area(compute_diameter(L, D, d, L))));
            if (i == 1 || stress_diff < stress_diff_at_location)
            {
                applied_force_location = i;
                stress_diff_at_location = stress_diff;
            }

            spring_system.forces(i) = 0; // reset force
        }

        // Progress bar update
        if (i % (segments / progress_bar_width) == 0)
        {
            std::cout << "#";
            std::cout.flush();
        }
    }
    std::cout << "] Done.\n";

    // now we print the distance where the force should be applied and the stress difference at that location
    std::cout << "Optimal force application location: Node " << applied_force_location
              << " (Position: " << nodes[static_cast<int>(applied_force_location)].position << " meters)" << std::endl;
    std::cout << "Minimum stress difference at this location: " << stress_diff_at_location << " Pa" << std::endl;

    sf::RenderWindow window(sf::VideoMode(window_size), "Spring System Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    sf::Clock deltaClock;

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

        ImGui::SFML::Update(window, deltaClock.restart());

        ImGui::Begin("System Controls");

        // add sliders to control scale_x and scale_y
        ImGui::SliderFloat("Scale X", &scale_x, 50.0f, 1000.0f);
        ImGui::SliderFloat("Scale Y", &scale_y, 1.0f, 10.0f);

        ImGui::Separator();

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                    1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        // add frametime graph
        static float frametimes[100] = {0};
        static int frametime_index = 0;
        frametimes[frametime_index] = 1000.0f / ImGui::GetIO().Framerate;
        frametime_index = (frametime_index + 1) % 100;
        ImGui::PlotLines("Frame Time (ms)", frametimes, 100, frametime_index, NULL, 0.0f, 50.0f, ImVec2(0, 80));

        ImGui::End();

        window.clear(sf::Color::Black);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}