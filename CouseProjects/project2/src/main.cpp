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

// Include project headers
#include "node.h"
#include "spring.h"
#include "spring_system.h"
#include "graphics.h"

float worldWidth = 1.2f;  // in meters
float worldHeight = 1.2f; // in meters

sf::Vector2f viewCenter(worldWidth - 0.5, worldHeight / 2.0);
float panSpeed = 0.01f;
float zoom = 1.0f;
sf::Vector2i lastMousePos;
bool isDragging = false;

int main()
{
    sf::Vector2u window_size(1200, 600);

    sf::RenderWindow window(sf::VideoMode(window_size), "2D Spring System Visualization");
    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    // Create a view for high-DPI / automatic scaling
    sf::View view;
    view.setSize(sf::Vector2f(worldWidth, -worldHeight));
    view.setCenter(sf::Vector2f(worldWidth / 2.f, -worldHeight / 2.f));
    window.setView(view);

    GraphicsRenderer renderer;

    sf::Clock deltaClock;

    const double A = 6.0e-4; // Cross-sectional area in m^2
    const double E = 210e9;  // Young's modulus in Pa

    // Define springs
    std::vector<Spring> springs = {
        Spring(0, 0, 1, A, E),
        Spring(1, 1, 2, A, E),
        Spring(2, 2, 0, A, E),
    };

    // Define nodes
    std::vector<Node> nodes = {
        Node(0, 0.0f, 0.0f, Fixed),           // Fixed node at origin
        Node(1, 0.0f, 1.0f, Slider, 90.0f),   // Slider node - vertical slider (allows up/down movement)
        Node(2, 0.577f, 1.0f, Slider, 60.0f), // Slider node - 60° slider (allows movement along 60° line)
    };

    // std::vector<Node> nodes = {
    //     Node(0, 0.0f, 0.0f, Fixed),
    //     Node(1, 1.0f, 0.0f, Free),
    //     Node(2, 1.0f, 0.57735f, Free),
    //     Node(3, 0.0f, 0.57735f, Fixed),
    // };

    // const double A = 6.0e-4;
    // const double E = 210e9;

    // std::vector<Spring> springs = {
    //     Spring(0, 0, 1, A, E),
    //     Spring(1, 1, 2, A, E),
    //     Spring(2, 2, 0, A, E),
    //     Spring(3, 2, 3, A, E),
    // };

    std::cout << "\n=== NODE CONFIGURATION ===" << std::endl;
    for (const auto &node : nodes)
    {
        std::cout << "Node " << node.id << ": position=(" << node.position[0] << ", " << node.position[1] << "), ";
        if (node.constraint_type == Fixed)
            std::cout << "FIXED";
        else if (node.constraint_type == Slider)
            std::cout << "SLIDER at " << node.constraint_angle << "°";
        else
            std::cout << "FREE";
        std::cout << std::endl;
    }

    SpringSystem spring_system(nodes, springs);

    // Apply a 400 kN downward force at node 2
    spring_system.forces(2 * 2 + 1) = -400000.0f; // Fy at node 2 (negative = downward)
    // spring_system.forces(2 * 2) = 400000.0f;

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

            // Mouse wheel zooming
            if (!ImGui::GetIO().WantCaptureMouse)
            {
                if (const auto *wheel = event->getIf<sf::Event::MouseWheelScrolled>())
                {
                    std::cout << "Mouse wheel delta: " << wheel->delta << std::endl;
                    if (wheel->delta > 0)
                        zoom = std::max(0.01f, zoom * 0.9f); // zoom in (clamped)
                    else
                        zoom = std::min(100.0f, zoom * 1.1f); // zoom out (clamped)
                }
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        ImGui::Begin("System Controls");

        ImGui::Text("2D Truss System Solver");
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
            sf::Color color = renderer.getStressColor(spring.stress, spring_system.min_stress, spring_system.max_stress);
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
        renderer.drawSystem(window, spring_system);

        // Mouse-based panning
        if (!ImGui::GetIO().WantCaptureMouse)
        {
            if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
            {
                sf::Vector2i mousePos = sf::Mouse::getPosition(window);

                if (!isDragging)
                {
                    lastMousePos = mousePos;
                    isDragging = true;
                }

                // Convert pixel delta into world delta using the current view
                sf::Vector2f worldLast = window.mapPixelToCoords(lastMousePos);
                sf::Vector2f worldNow = window.mapPixelToCoords(mousePos);

                sf::Vector2f delta = worldLast - worldNow;

                viewCenter += delta;
                lastMousePos = mousePos;
            }
            else
            {
                isDragging = false;
            }
        }

        // -------------------------
        // Aspect-correct view
        // -------------------------
        sf::Vector2u windowSize = window.getSize();
        float windowAspect = static_cast<float>(windowSize.x) / windowSize.y;
        float worldAspect = worldWidth / worldHeight;

        sf::Vector2f viewSize;
        if (windowAspect >= worldAspect)
        {
            viewSize.y = worldHeight * zoom;
            viewSize.x = viewSize.y * windowAspect;
        }
        else
        {
            viewSize.x = worldWidth * zoom;
            viewSize.y = viewSize.x / windowAspect;
        }

        sf::View view;
        view.setCenter(viewCenter);
        view.setSize(sf::Vector2f(viewSize.x, -viewSize.y)); // negative y = y-up
        window.setView(view);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}