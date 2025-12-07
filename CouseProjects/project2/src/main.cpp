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

    const double PI = 3.14159265358979323846;

    const double diameter_steel = 0.5;    // in
    const double diameter_aluminum = 0.4; // in

    const double A_steel = PI * std::pow(diameter_steel / 2.0, 2);          // in^2
    const double A_aluminum = PI * std::pow(diameter_aluminum / 2.0, 2);    // in^2
    const double I_steel = (PI / 64.0) * std::pow(diameter_steel, 4);       // in^4
    const double I_aluminum = (PI / 64.0) * std::pow(diameter_aluminum, 4); // in^4
    const double S_steel = I_steel / (diameter_steel / 2.0);                // in^3
    const double S_aluminum = I_aluminum / (diameter_aluminum / 2.0);       // in^3
    const double E_steel = 30e6;                                            // Psi
    const double E_aluminum = 10e6;                                         // Psi

    // idiot check print out material and beam properties
    std::cout << "Steel Material Properties:" << std::endl;
    std::cout << "  Young's Modulus: " << E_steel << " Psi" << std::endl;
    std::cout << "  Area: " << A_steel << " in^2" << std::endl;
    std::cout << "  Moment of Inertia: " << I_steel << " in^4" << std::endl;

    std::cout << "Aluminum Material Properties:" << std::endl;
    std::cout << "  Young's Modulus: " << E_aluminum << " Psi" << std::endl;
    std::cout << "  Area: " << A_aluminum << " in^2" << std::endl;
    std::cout << "  Moment of Inertia: " << I_aluminum << " in^4" << std::endl;

    MaterialProfile steel_material = {"Steel", E_steel};
    MaterialProfile aluminum_material = {"Aluminum", E_aluminum};
    BeamProfile steel_beam = {"Steel Beam", A_steel, I_steel, S_steel};
    BeamProfile aluminum_beam = {"Aluminum Beam", A_aluminum, I_aluminum, S_aluminum};
    BeamProfile steel_truss = {"Steel Truss", A_steel / 2.0f, 0.0f, 0.0f};
    BeamProfile aluminum_truss = {"Aluminum Truss", A_aluminum / 2.0f, 0.0f, 0.0f};

    // Define material profiles
    std::vector<MaterialProfile> material_profiles = {
        steel_material,
        aluminum_material,
    };

    // Define beam properties
    std::vector<BeamProfile> beam_profiles = {
        steel_beam,
        aluminum_beam,
        steel_truss,
        aluminum_truss,
    };

    std::vector<Node> nodes = {
        Node(0.0f, 0.0f, Fixed),
        Node(10.0f, 0.0f, Free),
    };

    std::vector<Beam> beams = {
        Beam(0, 1, 1, 1),
    };

    // Define nodes
    // std::vector<Node> nodes = {
    //     Node(12.0f, 0.0f, Free),
    //     Node(12.0f, 6.0f, Free),
    //     Node(0.0f, 0.0f, Slider, 90.0f),
    //     Node(0.0f, 10.0f, Free),
    //     Node(-10.0f, 10.0f, Fixed),
    // };

    // // Define beams (beams)
    // std::vector<Beam>
    //     beams = {
    //         Beam(0, 1, 0, 2),
    //         Beam(0, 2, 1, 3),
    //         Beam(1, 2, 0, 2),
    //         Beam(1, 3, 1, 3),
    //         Beam(2, 3, 0, 2),
    //         Beam(3, 4, 1, 1),
    //     };

    FEMSystem fem_system(nodes, beams, material_profiles, beam_profiles);

    fem_system.unit_system = Imperial;
    fem_system.forces(1 * 3 + 1) = -10.0; // Apply a downward force of 10 lbs at node 0

    // Course Project 2 Loads
    // fem_system.forces(0 * 3) = -1000.0;
    // fem_system.forces(0 * 3 + 1) = -1732.0;

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