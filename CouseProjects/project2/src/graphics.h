#pragma once
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include "spring_system.h"

class GraphicsRenderer
{
private:
    // Global dimensions for the world
    float worldWidth;
    float worldHeight;

    // View management
    sf::Vector2f viewCenter;
    float zoom;

    // Panning state
    sf::Vector2i lastMousePos;
    bool isDragging;

    // Window Focus
    bool isFocused;
    bool dragStartedInside;

    // Font for text rendering
    sf::Font font;
    mutable sf::Text text;

    void
    drawThickLine(sf::RenderTarget &target,
                  const sf::Vector2f &a,
                  const sf::Vector2f &b,
                  float thickness,
                  const sf::Color &color) const;

    void drawGrid(sf::RenderWindow &window) const;
    float getViewScale(const sf::RenderWindow &window) const;

public:
    GraphicsRenderer();

    // Initialize world dimensions
    void initialize(float width, float height);

    // DPI scaling
    float GetDPIScale(sf::WindowHandle handle);

    // Event handling
    void handleEvent(sf::RenderWindow &window, const sf::Event &event, bool mouseCapture);

    // Update panning (call each frame)
    void updatePanning(sf::RenderWindow &window, bool mouseCapture);

    // Update and apply view (call before drawing)
    void updateView(sf::RenderWindow &window);

    // Draw the spring system
    void drawSystem(sf::RenderWindow &window, const SpringSystem &system) const;

    // Get stress color
    sf::Color getStressColor(float stress, float min_stress, float max_stress) const;
};
