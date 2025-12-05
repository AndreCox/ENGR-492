#pragma once
#include "fem_system.h"
#include "graphics.h"
#include <imgui.h>
#include <SFML/Graphics.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Window/Keyboard.hpp>

class GUIHandler
{
public:
    GUIHandler(FEMSystem &system, GraphicsRenderer &renderer, sf::RenderWindow &window);

    void processEvent(const sf::Event &event);
    void render();

private:
    FEMSystem &fem_system;
    GraphicsRenderer &renderer;
    sf::RenderWindow &window;
    sf::Clock deltaClock;

    void applyDPIScale(float scale);
    void reloadFontsForDPI(float dpiScale);
    float pending_dpi_scale = 0.0f; // 0.0f means no change requested

    void systemControls();
    void nodeEditor();
    void beamEditor();
    void headerBar();
    void handleSavePopup();
    void handleLoadPopup();
    void handleDPIAdjust();

    // material properties editors
    void materialEditor();
    void profileEditor();

    void helpPage();

    bool show_system_controls = true;
    bool show_node_editor = false;
    bool show_beam_editor = false;
    bool request_save_popup = false;
    bool request_load_popup = false;
    bool request_dpi_adjust = false;
    bool show_help_page = false;

    char filename_buf[512] = "system";
    bool trigger_save_write = false;
    bool trigger_load_read = false;
    bool save_error = false;
    bool load_error = false;
    std::string error_msg = "";

    // material properties boolean flags
    bool show_material_editor = false;
    bool show_profile_editor = false;

    // DPI scaling members
    ImGuiStyle base_imgui_style;
    float base_font_global_scale = 1.0f;
    float current_dpi_scale = 1.0f; // 1.0 == 100%
};