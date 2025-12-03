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
    // Used for high dpi scaling and centering
    float worldWidth;
    float worldHeight;

    void drawThickLine(sf::RenderTarget &target,
                       const sf::Vector2f &a,
                       const sf::Vector2f &b,
                       float thickness,
                       const sf::Color &color) const;

public:
    GraphicsRenderer();

    float GetDPIScale(sf::WindowHandle handle);

    void setView(sf::RenderWindow &window);

    void drawSystem(sf::RenderWindow &window, const SpringSystem &system) const;

    sf::Color getStressColor(float stress, float min_stress, float max_stress) const;
};
