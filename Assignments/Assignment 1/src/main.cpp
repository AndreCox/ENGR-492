#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include <thread>

double shaft_radius(double x)
{
    return 10.0 + ((50 - 10) / (2000.0)) * x;
}

double elongation(double E, double F, double L, double A)
{
    return (F * L) / (A * E);
}

float scale_x = 0.4;
float scale_y = 4;
float scales[2] = {scale_x, scale_y};

int num_sections = 20;
float estimated_elongation = 0.0;

int main()
{
    sf::Vector2u window_size(1200, 400);
    const double shaft_length = 2000.0; // mm

    const double E = 210e9;        // Pa
    const double F = 950 * 1e3;    // N
    const double L = shaft_length; // mm

    sf::RenderWindow window(sf::VideoMode(window_size), "Shaft Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    // Scaling factors
    double min_diameter = shaft_radius(0);
    double max_diameter = shaft_radius(shaft_length);

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

        // 2. Show a simple window with various widgets

        ImGui::Begin("Shaft Controls");

        ImGui::Text("Shaft Length: %.1f mm", shaft_length);
        ImGui::Text("Min Diameter: %.1f mm", min_diameter);
        ImGui::Text("Max Diameter: %.1f mm", max_diameter);
        ImGui::Text("Young's Modulus: %.2f GPa", E / 1e9);
        ImGui::Text("Applied Force: %.2f kN", F / 1e3);

        ImGui::SliderInt("Number of Sections", &num_sections, 1, 5000);
        // add sliders for scale_x and scale_y
        ImGui::SliderFloat2("Scale (X, Y)", scales, 0.1f, 4.0f);
        ImGui::Text("Estimated Elongation: %.10f mm", estimated_elongation * 1e3);
        estimated_elongation = 0.0;

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
        double x_center = window_size.x / 2.0;
        sf::VertexArray shaft(sf::PrimitiveType::TriangleFan, 4);

        shaft[0].position = sf::Vector2f(x_center - (shaft_length * scales[0]) / 2, y_center - (shaft_radius(0) * scales[1]));
        shaft[1].position = sf::Vector2f(x_center + (shaft_length * scales[0]) / 2, y_center - (shaft_radius(shaft_length) * scales[1]));
        shaft[2].position = sf::Vector2f(x_center + (shaft_length * scales[0]) / 2, y_center + (shaft_radius(shaft_length) * scales[1]));
        shaft[3].position = sf::Vector2f(x_center - (shaft_length * scales[0]) / 2, y_center + (shaft_radius(0) * scales[1]));

        shaft[0].color = sf::Color::Blue;
        shaft[1].color = sf::Color::Blue;
        shaft[2].color = sf::Color::Blue;
        shaft[3].color = sf::Color::Blue;

        window.draw(shaft);

        // Draw sections of the shaft
        std::vector<sf::RectangleShape> sections;
        for (int i = 0; i < num_sections; ++i)
        {
            double x = (shaft_length / num_sections) * i;
            double radius = shaft_radius(x);
            sf::RectangleShape section(sf::Vector2f((shaft_length / num_sections) * scales[0], radius * 2 * scales[1]));
            section.setPosition(sf::Vector2f(
                static_cast<float>(x_center - (shaft_length * scales[0]) / 2 + (shaft_length / num_sections) * i * scales[0]),
                static_cast<float>(y_center - radius * scales[1])));

            section.setFillColor(sf::Color(250, 0, 0, 150));
            if (num_sections < 150)
            {
                section.setOutlineThickness(2);            // Set outline thickness
                section.setOutlineColor(sf::Color::Black); // Set outline color
            }
            sections.push_back(section);
        }

        for (int i = 0; i < num_sections; ++i)
        {
            double x = (shaft_length / num_sections) * i;
            double next_x = (shaft_length / num_sections) * (i + 1);
            double radius = shaft_radius(x);
            double next_radius = shaft_radius(next_x);
            double avg_radius = 0.5 * (radius + next_radius);
            double area = M_PI * avg_radius * avg_radius * 1e-6;          // mm^2 to m^2
            double segment_length = (shaft_length / num_sections) * 1e-3; // mm to m
            estimated_elongation += elongation(E, F, segment_length, area);
        }

        for (const auto &section : sections)
        {
            window.draw(section);
        }

        // allow user to scroll in and out with a 100ms buffer
        static sf::Clock keyBufferClock;
        if (keyBufferClock.getElapsedTime().asMilliseconds() > 100)
        {
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Up))
            {
                num_sections += 1;
                keyBufferClock.restart();
            }
            else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Down))
            {
                num_sections -= 1;
                if (num_sections < 1)
                    num_sections = 1;
                keyBufferClock.restart();
            }
        }

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
