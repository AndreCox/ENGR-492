#include "graphics.h"

void GraphicsRenderer::setView(sf::RenderWindow &window)
{
    // Create a view that shows the whole world rectangle
    sf::View view;
    view.setSize(sf::Vector2f(worldWidth, worldHeight));
    view.setCenter(sf::Vector2f(worldWidth / 2.0f, worldHeight / 2.0f));
    window.setView(view);
}

sf::Color GraphicsRenderer::getStressColor(float stress, float min_stress, float max_stress) const
{
    float normalized;
    float abs_max = std::max(std::abs(min_stress), std::abs(max_stress));

    if (abs_max < 1e-6f)
    {
        return sf::Color::White;
    }

    normalized = stress / abs_max;

    if (normalized < 0)
    {
        float t = -normalized;
        uint8_t r = static_cast<uint8_t>(255 * (1 - t));
        uint8_t g = static_cast<uint8_t>(255 * (1 - t));
        uint8_t b = 255;
        return sf::Color(r, g, b);
    }
    else
    {
        float t = normalized;
        uint8_t r = 255;
        uint8_t g = static_cast<uint8_t>(255 * (1 - t));
        uint8_t b = static_cast<uint8_t>(255 * (1 - t));
        return sf::Color(r, g, b);
    }
}

GraphicsRenderer::GraphicsRenderer()
{
}

void GraphicsRenderer::drawThickLine(sf::RenderTarget &target,
                                     const sf::Vector2f &a,
                                     const sf::Vector2f &b,
                                     float thickness,
                                     const sf::Color &color) const
{
    sf::Vector2f dir = b - a;
    float length = std::sqrt(dir.x * dir.x + dir.y * dir.y);
    if (length < 1e-6f)
        return;

    dir /= length;
    sf::Vector2f normal(-dir.y, dir.x);
    sf::Vector2f offset = normal * (thickness * 0.5f);

    sf::VertexArray quad(sf::PrimitiveType::Triangles, 6);

    quad[0].position = a + offset;
    quad[1].position = b + offset;
    quad[2].position = b - offset;

    quad[3].position = a + offset;
    quad[4].position = b - offset;
    quad[5].position = a - offset;

    for (int i = 0; i < 6; i++)
        quad[i].color = color;

    target.draw(quad);
}

void GraphicsRenderer::drawSystem(sf::RenderWindow &window, const SpringSystem &system) const
{
    // -------------------------
    // Draw springs with stress-based colors
    // -------------------------
    for (const auto &spring : system.springs)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        sf::Color springColor = getStressColor(spring.stress, system.min_stress, system.max_stress);

        sf::Vector2f p1(
            system.nodes[n1].position[0] + system.displacement(n1 * 2),
            system.nodes[n1].position[1] + system.displacement(n1 * 2 + 1));

        sf::Vector2f p2(
            system.nodes[n2].position[0] + system.displacement(n2 * 2),
            system.nodes[n2].position[1] + system.displacement(n2 * 2 + 1));

        // thickness in world units
        float thickness = 0.01f;

        drawThickLine(window, p1, p2, thickness, springColor);
    }

    // -------------------------
    // Draw nodes
    // -------------------------
    for (const auto &node : system.nodes)
    {
        sf::Vector2f pos(node.position[0] + system.displacement(node.id * 2),
                         node.position[1] + system.displacement(node.id * 2 + 1));

        if (node.constraint_type == Fixed)
        {
            float size = 0.025f;                                 // world units
            sf::RectangleShape square(sf::Vector2f(size, size)); // world units
            square.setOrigin(sf::Vector2f(size / 2, size / 2));  // center the square
            square.setPosition(pos);
            square.setFillColor(sf::Color::Red);
            window.draw(square);
        }
        else if (node.constraint_type == Slider)
        {
            sf::ConvexShape triangle(3);

            float perp_angle_rad = (node.constraint_angle + 90.0f) * M_PI / 180.0f;
            float dx = std::cos(perp_angle_rad);
            float dy = std::sin(perp_angle_rad);

            float size = 0.025f;                      // world units
            triangle.setPoint(0, sf::Vector2f(0, 0)); // tip
            triangle.setPoint(1, sf::Vector2f(-dx * size - dy * size / 2, -dy * size + dx * size / 2));
            triangle.setPoint(2, sf::Vector2f(-dx * size + dy * size / 2, -dy * size - dx * size / 2));

            triangle.setPosition(pos);
            triangle.setFillColor(sf::Color::Yellow);
            window.draw(triangle);
        }
        else
        {
            sf::CircleShape circle(0.025f); // radius in world units
            circle.setOrigin(sf::Vector2f(0.5f, 0.5f));
            circle.setPosition(pos);
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);
        }
    }

    // -------------------------
    // Draw forces as arrows
    // -------------------------
    for (int i = 0; i < system.nodes.size(); ++i)
    {
        const Node &node = system.nodes[i];

        // scale factor for force visualization (increase to make arrows shorter)
        const float arrowScale = 1000000.0f;
        float fx = system.forces(i * 2) / arrowScale;
        float fy = system.forces(i * 2 + 1) / arrowScale;

        if (std::abs(fx) < 1e-6f && std::abs(fy) < 1e-6f)
            continue;

        sf::Vector2f start(node.position[0] + system.displacement(i * 2),
                           node.position[1] + system.displacement(i * 2 + 1));
        sf::Vector2f end = start + sf::Vector2f(fx, fy);

        sf::Vertex line[2];
        line[0].position = start;
        line[0].color = sf::Color::Magenta;
        line[1].position = end;
        line[1].color = sf::Color::Magenta;
        window.draw(line, 2, sf::PrimitiveType::Lines);

        // Arrowhead
        sf::Vector2f dir = end - start;
        float length = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        if (length > 0)
        {
            sf::Vector2f unit = dir / length;
            sf::Vector2f perp(-unit.y, unit.x);
            float arrow_size = 0.025f; // world units

            sf::Vertex arrow[3];
            arrow[0].position = end;
            arrow[0].color = sf::Color::Magenta;
            arrow[1].position = end - unit * arrow_size + perp * (arrow_size / 2);
            arrow[1].color = sf::Color::Magenta;
            arrow[2].position = end - unit * arrow_size - perp * (arrow_size / 2);
            arrow[2].color = sf::Color::Magenta;

            window.draw(arrow, 3, sf::PrimitiveType::Triangles);
        }
    }
}