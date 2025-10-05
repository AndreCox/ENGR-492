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

double compute_springyness(double G, double J, double L)
{
    return (G * J) / L;
}

double compute_J(double diameter)
{
    return (M_PI * std::pow(diameter, 4)) / 32.0;
}

double compute_shear_stress(double torque, double radius, double J)
{
    return (torque * radius) / J;
}

class SpringSystem
{
public:
    std::vector<Node> nodes;
    std::vector<Spring> springs;
    Eigen::MatrixXf global_k_matrix;
    Eigen::VectorXf forces;
    Eigen::VectorXf angles; // angles in radians

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

        // Solve for angles (in radians)
        Eigen::VectorXf reduced_angles = reduced_k_matrix.colPivHouseholderQr().solve(reduced_forces);
        for (int i = 0; i < num_free_nodes; ++i)
        {
            angles(free_node_indices[i]) = reduced_angles(i);
        }

        // solve for forces
        forces = global_k_matrix * angles;

        std::cout << "Angles (radians):\n"
                  << angles << std::endl;

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
        angles = Eigen::VectorXf::Zero(num_nodes);

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

    const double G_AB = 3.4e6;      // lb/in^2
    const double G_BC = 4.0e6;      // lb/in^2
    const double D = 2.0;           // in
    const double d = 1.2;           // in
    const double JD = compute_J(D); // in^4
    const double Jd = compute_J(d); // in^4
    double T1 = -200;               // lb*ft
    double T2 = -100;               // lb*ft
    T1 *= 12;                       // lb*in
    T2 *= 12;                       // lb*in

    const double AD = 2.5 * 12; // in
    const double DB = 1.5 * 12; // in
    const double AB = AD + DB;  // in
    const double BC = 4.0 * 12; // in

    // Define nodes
    std::vector<Node> nodes = {
        Node(0, 0.0f, true), // Fixed node at position 0
        Node(1, AD, false),  // Free node at position AB
        Node(2, DB, false),  // Free node at position AB
        Node(3, BC, false),  // Free node at position AB
    };

    // Define springs
    std::vector<Spring> springs = {
        Spring(0, 0, 1, compute_springyness(G_AB, JD, AD)), // Spring between node 0 and 1
        Spring(1, 1, 2, compute_springyness(G_AB, JD, DB)), // Spring between node 1 and 2
        Spring(2, 2, 3, compute_springyness(G_BC, Jd, BC)), // Spring between node 2 and 3
    };

    SpringSystem spring_system(nodes, springs);
    spring_system.forces << 0, T1, 0, T2; // Apply torques at nodes 1 and 3
    spring_system.solve_system();

    // determine the shear stress at node 0 and node 2
    double shear_stress_node_0 = compute_shear_stress(spring_system.forces(0), D / 2.0, JD);
    double torque_in_spring_2 = springs[2].k * (spring_system.angles(3) - spring_system.angles(2));
    double shear_stress_node_2 = compute_shear_stress(std::abs(torque_in_spring_2), d / 2.0, Jd);
    std::cout << "Shear Stress at Node 0: " << shear_stress_node_0 << " lb/in^2" << std::endl;
    std::cout << "Shear Stress at Node 2: " << shear_stress_node_2 << " lb/in^2" << std::endl;
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