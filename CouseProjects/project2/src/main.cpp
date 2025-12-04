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

int main()
{
    sf::Vector2u window_size(1200, 600);
    sf::RenderWindow window(sf::VideoMode(window_size), "2D Spring System Visualization");

    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    GraphicsRenderer renderer;

    // Initialize the renderer
    renderer.initialize(1.2f, 1.2f);

    ImGuiIO &io = ImGui::GetIO();
    float dpiScale = renderer.GetDPIScale(window.getNativeHandle());
    io.FontGlobalScale = 1.8f;             // Scale text
    ImGui::GetStyle().ScaleAllSizes(1.8f); // Scale widgets

    sf::Clock deltaClock;

    // Course Project 2 Example
    std::vector<Node> nodes = {
        Node(12.0f, 0.0f, Free),
        Node(12.0f, 6.0f, Free),
        Node(0.0f, 0.0f, Slider, 90.0f),
        Node(0.0f, 10.0f, Fixed),
    };

    const double A_steel = 0.1963;    // in^2
    const double A_aluminum = 0.1257; // in^2
    const double E_steel = 30e6;      // Psi
    const double E_aluminum = 10e6;   // Psi

    std::vector<Spring> springs = {
        Spring(0, 1, A_steel, E_steel),
        Spring(0, 2, A_aluminum, E_aluminum),
        Spring(1, 2, A_steel, E_steel),
        Spring(1, 3, A_aluminum, E_aluminum),
        Spring(2, 3, A_steel, E_steel),
    };

    std::cout << "\n=== NODE CONFIGURATION ===" << std::endl;
    int index = 0;
    for (const auto &node : nodes)
    {
        std::cout << "Node " << index << ": position=(" << node.position[0] << ", " << node.position[1] << "), ";
        if (node.constraint_type == Fixed)
            std::cout << "FIXED";
        else if (node.constraint_type == Slider)
            std::cout << "SLIDER at " << node.constraint_angle << "°";
        else
            std::cout << "FREE";
        std::cout << std::endl;
        ++index;
    }

    SpringSystem spring_system(nodes, springs);

    // Course Project 2 Loads
    spring_system.forces(0 * 2) = -1000.0;
    spring_system.forces(0 * 2 + 1) = -1732.0;

    // Apply a 400 kN downward force at node 2
    // spring_system.forces(2 * 2 + 1) = -4.0e5; // Fy at node 1 (negative = downward)
    // spring_system.forces(1 * 2 + 1) = -5.0e5; // Fy at node 2 (negative = downward)
    // spring_system.forces(2 * 2) = 4.0e5;      // Fy at node 3 (negative = downward)

    std::cout
        << "\nApplied Forces (N):" << std::endl;
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
            renderer.handleEvent(window, *event, ImGui::GetIO().WantCaptureMouse);
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
        int index = 0;
        for (const auto &spring : spring_system.springs)
        {
            sf::Color color = renderer.getStressColor(spring.stress, spring_system.min_stress, spring_system.max_stress);
            ImGui::TextColored(ImVec4(color.r / 255.0f, color.g / 255.0f, color.b / 255.0f, 1.0f),
                               "Spring %d: %.2f MPa", index, spring.stress);
            ++index;
        }

        ImGui::Separator();
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                    1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        ImGui::End();

        window.clear(sf::Color(30, 30, 30));

        // Update panning
        renderer.updatePanning(window, ImGui::GetIO().WantCaptureMouse);

        // Update and apply view
        renderer.updateView(window);

        // Draw the spring system
        renderer.drawSystem(window, spring_system);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}