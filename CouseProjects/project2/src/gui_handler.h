#pragma once
#include "spring_system.h"
#include "graphics.h"
#include <imgui.h>
#include <SFML/Graphics.hpp>
#include <SFML/Window/Event.hpp>

class GUIHandler
{
public:
    GUIHandler(SpringSystem &system, GraphicsRenderer &renderer, sf::RenderWindow &window);

    void processEvent(const sf::Event &event);
    void render();

private:
    SpringSystem &spring_system;
    GraphicsRenderer &renderer;
    sf::RenderWindow &window;
    sf::Clock deltaClock;

    void systemControls();
    void nodeEditor();
    void springEditor();
    void headerBar();
    void handleSavePopup();
    void handleLoadPopup();

    bool show_system_controls = true;
    bool show_node_editor = false;
    bool show_spring_editor = false;
    bool request_save_popup = false;
    bool request_load_popup = false;
};