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
#include "beam.h"
#include "fem_system.h"
#include "graphics.h"
#include "gui_handler.h"

int main()
{
    sf::Vector2u window_size(1200, 600);
    sf::RenderWindow window(sf::VideoMode(window_size), "Fast FEM");

    if (!ImGui::SFML::Init(window))
    {
        std::cerr << "Failed to initialize ImGui-SFML" << std::endl;
        return -1;
    }

    const double A_steel = 0.1963;    // in^2
    const double A_aluminum = 0.1257; // in^2
    const double E_steel = 30e6;      // Psi
    const double E_aluminum = 10e6;   // Psi

    MaterialProfile steel_material = {"Steel", E_steel};
    MaterialProfile aluminum_material = {"Aluminum", E_aluminum};
    BeamProfile steel_profile = {"Steel Beam", A_steel};
    BeamProfile aluminum_profile = {"Aluminum Beam", A_aluminum};

    // Define material profiles
    std::vector<MaterialProfile> material_profiles = {
        steel_material,
        aluminum_material,
    };

    // Define beam properties
    std::vector<BeamProfile> beam_profiles = {
        steel_profile,
        aluminum_profile,
    };

    // Define nodes
    std::vector<Node> nodes = {
        Node(12.0f, 0.0f, Free),
        Node(12.0f, 6.0f, Free),
        Node(0.0f, 0.0f, Slider, 90.0f),
        Node(0.0f, 10.0f, FixedPin),
    };

    // Define beams (beams)
    std::vector<Beam>
        beams = {
            Beam(0, 1, 0, 0),
            Beam(0, 2, 1, 1),
            Beam(1, 2, 0, 0),
            Beam(1, 3, 1, 1),
            Beam(2, 3, 0, 0),
        };

    FEMSystem fem_system(nodes, beams, material_profiles, beam_profiles);

    // Course Project 2 Loads
    fem_system.forces(0 * 3) = -1000.0;
    fem_system.forces(0 * 3 + 1) = -1732.0;

    // Initialize the renderer
    GraphicsRenderer renderer(fem_system);
    renderer.initialize(1.1f, 1.1f);
    renderer.autoZoomToFit();

    // Create GUI Handler
    GUIHandler gui_handler(fem_system, renderer, window);

    fem_system.solve_system();

    while (window.isOpen())
    {
        while (const std::optional<sf::Event> event = window.pollEvent())
        {

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
            // Handle ImGui events
            gui_handler.processEvent(*event);
        }

        gui_handler.render();

        window.clear(sf::Color(30, 30, 30));

        // Update panning
        renderer.updatePanning(window, ImGui::GetIO().WantCaptureMouse);

        // Update and apply view
        renderer.updateView(window);

        // Draw the spring system
        renderer.drawSystem(window);

        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}