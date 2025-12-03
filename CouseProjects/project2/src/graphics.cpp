#include "graphics.h"

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

GraphicsRenderer::GraphicsRenderer(float scale_x, float scale_y)
    : scale_x(scale_x), scale_y(scale_y) {}

void GraphicsRenderer::setScale(float sx, float sy)
{
    scale_x = sx;
    scale_y = sy;
}

void GraphicsRenderer::drawSystem(sf::RenderWindow &window, const SpringSystem &system, float offset_x, float offset_y) const
{
    // Draw springs with stress-based colors
    for (const auto &spring : system.springs)
    {
        int n1 = spring.nodes[0];
        int n2 = spring.nodes[1];

        sf::Color springColor = getStressColor(spring.stress, system.min_stress, system.max_stress);

        sf::Vertex line[2];
        line[0].position = sf::Vector2f(offset_x + system.nodes[n1].position[0] * scale_x + system.displacement(n1 * 2) * scale_x,
                                        offset_y - system.nodes[n1].position[1] * scale_y - system.displacement(n1 * 2 + 1) * scale_y);
        line[0].color = springColor;

        line[1].position = sf::Vector2f(offset_x + system.nodes[n2].position[0] * scale_x + system.displacement(n2 * 2) * scale_x,
                                        offset_y - system.nodes[n2].position[1] * scale_y - system.displacement(n2 * 2 + 1) * scale_y);
        line[1].color = springColor;

        window.draw(line, 2, sf::PrimitiveType::Lines);
    }

    // Draw nodes
    for (const auto &node : system.nodes)
    {
        if (node.constraint_type == Fixed)
        {
            sf::RectangleShape square(sf::Vector2f(10, 10));
            square.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x - 5,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y - 5));
            square.setFillColor(sf::Color::Red);
            window.draw(square);
        }
        else if (node.constraint_type == Slider)
        {
            sf::ConvexShape triangle(3); // Create a triangle with 3 points
            triangle.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y));

            // Direction perpendicular to the slider (add 90° to the constraint angle)
            float perp_angle_rad = (node.constraint_angle + 90.0f) * M_PI / 180.0f;
            float dx = std::cos(perp_angle_rad);
            float dy = -std::sin(perp_angle_rad); // Negative because y-axis is inverted in SFML

            // Triangle pointing in the perpendicular direction
            triangle.setPoint(0, sf::Vector2f(0, 0));                                 // Tip of the triangle
            triangle.setPoint(1, sf::Vector2f(-dx * 10 - dy * 5, -dy * 10 + dx * 5)); // Base left
            triangle.setPoint(2, sf::Vector2f(-dx * 10 + dy * 5, -dy * 10 - dx * 5)); // Base right
            triangle.setFillColor(sf::Color::Yellow);
            window.draw(triangle);
        }
        else
        {
            sf::CircleShape circle(5);
            circle.setPosition(sf::Vector2f(
                offset_x + node.position[0] * scale_x + system.displacement(node.id * 2) * scale_x - 5,
                offset_y - node.position[1] * scale_y - system.displacement(node.id * 2 + 1) * scale_y - 5));
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);
        }
    }

    // Draw forces as arrows on nodes
    for (int i = 0; i < system.nodes.size(); ++i)
    {
        const Node &node = system.nodes[i];
        float fx = system.forces(i * 2) / 1000.0f;      // Scale down for visualization
        float fy = -system.forces(i * 2 + 1) / 1000.0f; // Invert y for SFML
        if (std::abs(fx) < 1e-3f && std::abs(fy) < 1e-3f)
            continue; // Skip zero forces
        sf::Vector2f start(offset_x + node.position[0] * scale_x + system.displacement(i * 2) * scale_x,
                           offset_y - node.position[1] * scale_y - system.displacement(i * 2 + 1) * scale_y);
        sf::Vector2f end = start + sf::Vector2f(fx, fy);
        sf::Vertex line[2];
        line[0].position = start;
        line[0].color = sf::Color::Magenta;
        line[1].position = end;
        line[1].color = sf::Color::Magenta;

        // Draw arrowhead
        sf::Vector2f direction = end - start;
        float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
        if (length > 0)
        {
            sf::Vector2f unit_dir = direction / length;
            sf::Vector2f perp_dir(-unit_dir.y, unit_dir.x);
            float arrow_size = 5.0f;
            sf::Vertex arrow[3];
            arrow[0].position = end;
            arrow[0].color = sf::Color::Magenta;
            arrow[1].position = end - unit_dir * arrow_size + perp_dir * (arrow_size / 2);
            arrow[1].color = sf::Color::Magenta;
            arrow[2].position = end - unit_dir * arrow_size - perp_dir * (arrow_size / 2);
            arrow[2].color = sf::Color::Magenta;
            window.draw(arrow, 3, sf::PrimitiveType::Triangles);
        }

        window.draw(line, 2, sf::PrimitiveType::Lines);
    }
}
