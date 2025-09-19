#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/Text.hpp>
#include <imgui.h>
#include <imgui-SFML.h>
#include <Eigen/Eigen>
#include <vector>

float scale_x = 1000;
float scale_y = 4;
float scales[2] = {scale_x, scale_y};

// Shaft properties
float steel_length = 0.3;    // meters
float aluminum_length = 0.2; // meters

float steel_modulus = 200e9;   // Pascals
float aluminum_modulus = 70e9; // Pascals

float circular_area = 1.13097e-4; // m^2 (for a diameter of 12 mm)

class Shaft
{
public:
    float length;
    float modulus;
    float circular_area;
    Eigen::Matrix2f k_matrix; // stiffness matrix

    Shaft(float len, float mod, float area)
        : length(len), modulus(mod), circular_area(area)
    {
        compute_stiffness();
    }

private:
    double k; // stiffness
    void compute_stiffness()
    {
        k = (modulus * circular_area) / length;
        k_matrix << k, -k,
            -k, k;
    }
};

class ShaftSystem
{
public:
    std::vector<Shaft> shafts;
    Eigen::MatrixXf global_k_matrix;
    Eigen::VectorXf forces;
    Eigen::VectorXf displacements;

    ShaftSystem(const std::vector<Shaft> &shaft_list)
        : shafts(shaft_list)
    {
        forces = Eigen::VectorXf::Zero(shafts.size() + 1);
        displacements = Eigen::VectorXf::Zero(shafts.size() + 1);
        assemble_global_stiffness();
    }

    void solve_system(const Eigen::VectorXf &f, const Eigen::VectorXf &d)
    {
        forces = f;
        displacements = d;

        int size = shafts.size() + 1;
        Eigen::MatrixXf k_reduced = global_k_matrix;
        Eigen::VectorXf f_reduced = forces;
        Eigen::VectorXf d_reduced = displacements;

        // Apply boundary conditions
        for (int i = 0; i < size; ++i)
        {
            if (!std::isnan(displacements(i)))
            {
                // Known displacement, modify the system
                k_reduced.row(i).setZero();
                k_reduced.col(i).setZero();
                k_reduced(i, i) = 1.0f;
                f_reduced(i) = displacements(i);
            }
            else if (!std::isnan(forces(i)))
            {
                // Known force, do nothing
            }
            else
            {
                // Neither force nor displacement is known, this is an error
                std::cerr << "Error: Node " << i << " has neither force nor displacement defined." << std::endl;
                return;
            }
        }

        // Solve the reduced system
        Eigen::VectorXf d_solution = k_reduced.colPivHouseholderQr().solve(f_reduced);

        // Update displacements and forces
        for (int i = 0; i < size; ++i)
        {
            if (std::isnan(displacements(i)))
            {
                displacements(i) = d_solution(i);
            }
        }

        forces = global_k_matrix * displacements;

        // Output results
        std::cout << "Displacements (m):\n"
                  << displacements << std::endl;
        std::cout << "Forces (N):\n"
                  << forces << std::endl;
    }

private:
    void
    assemble_global_stiffness()
    {
        int size = shafts.size() + 1;
        global_k_matrix = Eigen::MatrixXf::Zero(size, size);
        for (size_t i = 0; i < shafts.size(); ++i)
        {
            global_k_matrix.block<2, 2>(i, i) += shafts[i].k_matrix;
        }

        std::cout << "Global Stiffness Matrix (N/m):\n"
                  << global_k_matrix << std::endl;
    }
};

int main()
{

    sf::Font font;
    if (!font.openFromFile("assets/roboto.ttf"))
    {
        std::cerr << "Failed to load font\n";
        return -1;
    }

    sf::Vector2u window_size(1200, 400);

    sf::RenderWindow window(sf::VideoMode(window_size), "Shaft Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    sf::Clock deltaClock;

    // Generate our 2 shafts that make up the system
    Shaft aluminum_shaft(aluminum_length, aluminum_modulus, circular_area);
    Shaft steel_shaft(steel_length, steel_modulus, circular_area);

    std::vector<Shaft> shafts = {aluminum_shaft, steel_shaft};
    ShaftSystem shaft_system(shafts);
    // Define boundary conditions
    Eigen::VectorXf forces(3);
    Eigen::VectorXf displacements(3);
    forces << NAN, 0, 100e3;      // Apply a force of 100 kN at the end of the aluminum shaft
    displacements << 0, NAN, NAN; // Fix the left end of the steel shaft

    shaft_system.solve_system(forces, displacements);

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

        ImGui::Begin("Shaft Controls");

        // add sliders to control scale_x and scale_y
        ImGui::SliderFloat("Scale X", &scale_x, 50.0f, 1000.0f);
        ImGui::SliderFloat("Scale Y", &scale_y, 1.0f, 10.0f);

        // allow user to modify force applied
        static float applied_force = 100e3;
        if (ImGui::InputFloat("Applied Force (N)", &applied_force))
        {
            forces(2) = applied_force;
            shaft_system.solve_system(forces, displacements);
        }

        // allow user to modify the material properties with sliders
        static float steel_E = steel_modulus / 1e9;
        static float aluminum_E = aluminum_modulus / 1e9;
        if (ImGui::SliderFloat("Steel Modulus (GPa)", &steel_E, 100.0f, 300.0f))
        {
            shafts[1] = Shaft(steel_length, steel_E * 1e9, circular_area);
            shaft_system = ShaftSystem(shafts);
            shaft_system.solve_system(forces, displacements);
        }
        if (ImGui::SliderFloat("Aluminum Modulus (GPa)", &aluminum_E, 50.0f, 150.0f))
        {
            shafts[0] = Shaft(aluminum_length, aluminum_E * 1e9, circular_area);
            shaft_system = ShaftSystem(shafts);
            shaft_system.solve_system(forces, displacements);
        }
        // length sliders
        static float steel_L = steel_length;
        static float aluminum_L = aluminum_length;
        if (ImGui::SliderFloat("Steel Length (m)", &steel_L, 0.1f, 1.0f))
        {
            shafts[1] = Shaft(steel_L, steel_E * 1e9, circular_area);
            shaft_system = ShaftSystem(shafts);
            shaft_system.solve_system(forces, displacements);
        }
        if (ImGui::SliderFloat("Aluminum Length (m)", &aluminum_L, 0.1f, 1.0f))
        {
            shafts[0] = Shaft(aluminum_L, aluminum_E * 1e9, circular_area);
            shaft_system = ShaftSystem(shafts);
            shaft_system.solve_system(forces, displacements);
        }

        ImGui::Separator();

        // display solved displacements and forces
        ImGui::Text("Displacements (m):");
        for (int i = 0; i < shaft_system.displacements.size(); ++i)
        {
            ImGui::Text("Node %d: %.6f", i, shaft_system.displacements(i));
        }

        ImGui::Separator();
        ImGui::Text("Forces (N):");
        for (int i = 0; i < shaft_system.forces.size(); ++i)
        {
            ImGui::Text("Node %d: %.2f", i, shaft_system.forces(i));
        }

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                    1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        // add frametime graph
        static float frametimes[100] = {0};
        static int frametime_index = 0;
        frametimes[frametime_index] = 1000.0f / ImGui::
                                                    GetIO()
                                                        .Framerate;
        frametime_index = (frametime_index + 1) % 100;
        ImGui::PlotLines("Frame Time (ms)", frametimes, 100, frametime_index, NULL, 0.0f, 50.0f, ImVec2(0, 80));

        ImGui::End();

        window.clear(sf::Color::Black);

        // Draw shaft in center of window
        double y_center = window_size.y / 2.0;
        double x_start = 100.0; // left margin

        // Draw shafts with deformation
        for (size_t i = 0; i < shafts.size(); ++i)
        {
            // Get node positions before and after deformation
            double node0_x = x_start + (i == 0 ? 0 : shafts[i - 1].length * scale_x + (i - 1 >= 0 ? 0 : 0));
            double node1_x = node0_x + shafts[i].length * scale_x;

            double disp0 = shaft_system.displacements(i) * scale_x;
            double disp1 = shaft_system.displacements(i + 1) * scale_x;
            double node0_x_disp = node0_x + disp0;
            double node1_x_disp = node1_x + disp1;
            double node0_y = y_center;
            double node1_y = y_center;

            // Draw shaft as a rectangle
            sf::RectangleShape shaft_shape(sf::Vector2f(node1_x_disp - node0_x_disp, 20));
            shaft_shape.setPosition(sf::Vector2f(static_cast<float>(node0_x_disp), static_cast<float>(node0_y - 10)));
            shaft_shape.setFillColor(i == 0 ? sf::Color::Green : sf::Color::Blue);
            window.draw(shaft_shape);

            // Draw nodes as circles
            sf::CircleShape node_shape(10);
            node_shape.setOrigin(sf::Vector2f(10.f, 10.f));
            node_shape.setPosition(sf::Vector2f(static_cast<float>(node0_x), static_cast<float>(node0_y)));
            node_shape.setFillColor(sf::Color::Red);
            window.draw(node_shape);
            if (i == shafts.size() - 1)
            {
                node_shape.setPosition(sf::Vector2f(static_cast<float>(node1_x), static_cast<float>(node1_y)));
                window.draw(node_shape);
            }
        }

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
