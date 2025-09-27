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

    void draw(sf::RenderWindow &window, const std::vector<Node> &nodes) const
    {
        // Draw a zigzag (resistor-like) spring between the two nodes
        const sf::Vector2f p1(nodes[this->nodes[0]].position * (scale_x / 1000.0f), window.getSize().y / 2);
        const sf::Vector2f p2(nodes[this->nodes[1]].position * (scale_x / 1000.0f), window.getSize().y / 2);

        const int num_zigs = 8;        // Number of zigzag segments
        const float amplitude = 15.0f; // Height of zigzag

        // Calculate direction vector and total length
        sf::Vector2f dir = p2 - p1;
        float length = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        sf::Vector2f unit_dir = dir / length;

        // Make straight segments proportional to total length (but with minimum)
        float straight_len = std::max(15.0f, length * 0.15f); // At least 15 pixels or 15% of total length

        std::vector<sf::Vertex> vertices;

        // Start straight segment
        sf::Vector2f straight_start = p1;
        sf::Vector2f zig_start = p1 + unit_dir * straight_len;

        // End positions
        sf::Vector2f zig_end = p2 - unit_dir * straight_len;
        sf::Vector2f straight_end = p2;

        // Add start straight line
        sf::Vertex startVertex;
        startVertex.position = straight_start;
        startVertex.color = sf::Color::White;
        vertices.push_back(startVertex);

        sf::Vertex zigStartVertex;
        zigStartVertex.position = zig_start;
        zigStartVertex.color = sf::Color::White;
        vertices.push_back(zigStartVertex);

        // Zigzag segments - interpolate between zig_start and zig_end
        for (int i = 1; i <= num_zigs; ++i)
        {
            float t = static_cast<float>(i) / (num_zigs + 1);
            sf::Vector2f pos = zig_start + t * (zig_end - zig_start);

            // Alternate the direction of the zigzag
            float offset = (i % 2 == 0 ? -amplitude : amplitude);

            // Perpendicular direction
            sf::Vector2f perp(-unit_dir.y, unit_dir.x); // Perpendicular to direction

            sf::Vertex vertex;
            vertex.position = pos + perp * offset;
            vertex.color = sf::Color::White;
            vertices.push_back(vertex);
        }

        // Add end straight line
        sf::Vertex zigEndVertex;
        zigEndVertex.position = zig_end;
        zigEndVertex.color = sf::Color::White;
        vertices.push_back(zigEndVertex);

        sf::Vertex endVertex;
        endVertex.position = straight_end;
        endVertex.color = sf::Color::White;
        vertices.push_back(endVertex);

        window.draw(vertices.data(), vertices.size(), sf::PrimitiveType::LineStrip);
    }

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

    void draw(sf::RenderWindow &window)
    {
        for (const auto &spring : springs)
        {
            spring.draw(window, nodes);
        }
        for (const auto &node : nodes)
        {
            node.draw(window);
        }
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

        std::cout << "Displacements:\n"
                  << displacements << std::endl;

        std::cout << "Forces:\n"
                  << forces << std::endl;
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

        std::cout << "Global Stiffness Matrix:\n"
                  << global_k_matrix << std::endl;
    }
};

int main()
{
    sf::Vector2u window_size(1200, 400);

    sf::RenderWindow window(sf::VideoMode(window_size), "Spring System Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    sf::Clock deltaClock;

    const double E = 10.6e6;         // Psi
    const double thickness = 0.0125; // in

    // system 1
    // // average widths down the system
    // const double W1 = compute_width(2.5 / 2.0);        // in
    // const double W2 = (compute_width(5) / 2.0) - 0.25; // in
    // const double W3 = (compute_width(8.75));           // in

    // // lengths of each segment
    // const double L1 = 2.5; // in
    // const double L2 = 5.0; // in
    // const double L3 = 2.5; // in

    // // compute the springs
    // const double k1 = compute_springyness(E, W1 * thickness, L1);
    // const double k23 = compute_springyness(E, W2 * thickness, L2);
    // const double k4 = compute_springyness(E, W3 * thickness, L3);

    // // Define nodes
    // std::vector<Node>
    //     nodes = {
    //         Node(0, 0.0f, true),    // Fixed node at position 0
    //         Node(1, 200.0f, false), // Free node at position 200
    //         Node(2, 400.0f, false), // Free node at position 400
    //         Node(3, 600.0f, false)  // Free node at position 600
    //     };

    // // Define springs
    // std::vector<Spring> springs = {
    //     Spring(0, 0, 1, k1),
    //     Spring(1, 1, 2, k23),
    //     Spring(2, 1, 2, k23),
    //     Spring(3, 2, 3, k4),
    // };

    // system 2
    // average widths down the system
    const double W1 = compute_width(2.5 / 2.0);
    const double W2 = (compute_width(2.5 + (2.5 / 2.0)) / 2.0) - 0.25;
    const double W3 = (compute_width(5 + (2.5 / 2.0)) / 2.0) - 0.25;
    const double W4 = compute_width(8.75);

    // lengths of each segment
    const double L1 = 2.5;
    const double L2 = 2.5;
    const double L3 = 2.5;
    const double L4 = 2.5;

    // compute the springs
    const double k1 = compute_springyness(E, W1 * thickness, L1);
    const double k24 = compute_springyness(E, W2 * thickness, L2);
    const double k35 = compute_springyness(E, W3 * thickness, L3);
    const double k6 = compute_springyness(E, W4 * thickness, L4);

    // Define nodes
    std::vector<Node> nodes = {
        Node(0, 0.0f, true),    // Fixed node at position 0
        Node(1, 200.0f, false), // Free node at position 200
        Node(2, 400.0f, false), // Free node at position 400
        Node(3, 400.0f, false), // Free node at position 400
        Node(4, 600.0f, false), // Free node at position 800
        Node(5, 800.0f, false)  // Free node at position 1000
    };

    // Define springs
    std::vector<Spring> springs = {
        Spring(0, 0, 1, k1),
        Spring(1, 1, 2, k24),
        Spring(2, 1, 3, k24),
        Spring(3, 2, 4, k35),
        Spring(4, 3, 4, k35),
        Spring(5, 4, 5, k6),
    };

    SpringSystem spring_system(nodes, springs);
    // spring_system.forces << 0, 0, 0, 1500; // Apply forces to free nodes
    spring_system.forces << 0, 0, 0, 0, 0, 1500; // Apply forces to free nodes
    spring_system.solve_system();

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

        // draw the spring system
        spring_system.draw(window);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}