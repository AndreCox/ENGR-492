#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include "spring_system.h"

class GraphicsRenderer
{
private:
    float scale_x;
    float scale_y;

public:
    GraphicsRenderer(float scale_x, float scale_y);

    void setScale(float sx, float sy);

    void drawSystem(sf::RenderWindow &window, const SpringSystem &system, float offset_x, float offset_y) const;

    sf::Color getStressColor(float stress, float min_stress, float max_stress) const;
};
