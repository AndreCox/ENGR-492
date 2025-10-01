#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp>
#include <imgui.h>
#include <imgui-SFML.h>
#include <Eigen/Eigen>
#include <vector>

float scale_x = 1000;
float scale_y = 4;
float scales[2] = {scale_x, scale_y};

float D = 0.75;  // major diameter in meters
float d = 0.25;  // minor diameter in meters
float L = 1;     // length in meters
float E = 260e9; // Young's modulus in Pascals
float P = 500;   // applied torque in Newton-meters

int num_segments = 10;

float calculate_area(float diameter)
{
    return M_PI * pow(diameter / 2.0, 2);
}

float calculate_stretch(float diameter, float length, float force, float youngs_modulus)
{
    float area = calculate_area(diameter);
    return (force * length) / (area * youngs_modulus);
}

float calculate_stress(float diameter, float force)
{
    float area = calculate_area(diameter);
    return force / area;
}

float calculate_spring(float diameter, float length, float youngs_modulus)
{
    float area = calculate_area(diameter);
    return (area * youngs_modulus) / length;
}

int main()
{

    // create a matrix to hold the spring global spring matrix should be the size of the number of segments times 2 (for x and y) by the number of segments times 2
    Eigen::MatrixXd global_stiffness_matrix = Eigen::MatrixXd::Zero(num_segments * 2, num_segments * 2);

    for (int i = 0; i < num_segments; ++i)
    {
        float segment_length = L / num_segments;
        float stretch = calculate_stretch(D, segment_length, P, E);
        float stress = calculate_stress(D, P);
        float spring_constant = calculate_spring(D, segment_length, E);

        // local stiffness matrix for the segment
        Eigen::Matrix2d local_stiffness_matrix;
        local_stiffness_matrix << spring_constant, -spring_constant,
            -spring_constant, spring_constant;

        // assemble into global stiffness matrix
        global_stiffness_matrix.block<2, 2>(i * 2, i * 2) += local_stiffness_matrix;
        if (i > 0)
        {
            global_stiffness_matrix.block<2, 2>(i * 2, (i - 1) * 2) += -local_stiffness_matrix;
            global_stiffness_matrix.block<2, 2>((i - 1) * 2, i * 2) += -local_stiffness_matrix;
        }
    }

    std::cout << "Global Stiffness Matrix:\n"
              << global_stiffness_matrix << std::endl;

    sf::Vector2u window_size(1200, 400);

    sf::RenderWindow window(sf::VideoMode(window_size), "Shaft Visualization");
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

        ImGui::Begin("Shaft Controls");

        // add sliders to control scale_x and scale_y
        ImGui::SliderFloat("Scale X", &scale_x, 50.0f, 1000.0f);
        ImGui::SliderFloat("Scale Y", &scale_y, 1.0f, 10.0f);

        ImGui::Separator();

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

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
