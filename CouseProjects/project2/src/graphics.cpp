#include "graphics.h"

GraphicsRenderer::GraphicsRenderer()
    : worldWidth(1.2f),
      worldHeight(1.2f),
      viewCenter(worldWidth / 2.0f, worldHeight / 2.0f),
      zoom(20.0f),
      isDragging(false),
      isFocused(true),
      dragStartedInside(false),
      font(),
      text(font)
{
    if (!font.openFromFile("resources/fonts/Roboto-Regular.ttf"))
    {
        std::cerr << "Failed to load font!" << std::endl;
    }

    text.setFont(font);

    text.setFillColor(sf::Color::Black);
    text.setOutlineColor(sf::Color::White); // outline color
    text.setOutlineThickness(4.0f);         // thickness in pixels
    text.setStyle(sf::Text::Bold);
    text.setCharacterSize(55); // initial size; can adjust per zoom
}

float GraphicsRenderer::GetDPIScale(sf::WindowHandle handle)
{
    // -----------------------------
    // WINDOWSW
    // -----------------------------
#if defined(_WIN32)
    UINT dpi = GetDpiForWindow(handle);
    return dpi / 96.0f;

    // -----------------------------
    // MACOS
    // -----------------------------
#elif defined(__APPLE__)
    // macOS always uses 72 DPI as base, scale is backingScaleFactor
    float scale = [[NSScreen mainScreen] backingScaleFactor];
    return scale;

    // -----------------------------
    // LINUX (X11)
    // -----------------------------
#elif defined(__linux__) && defined(SFML_SYSTEM_X11)
    Display *dpy = XOpenDisplay(nullptr);
    if (!dpy)
        return 1.0f;

    char *res = XResourceManagerString(dpy);
    if (!res)
    {
        XCloseDisplay(dpy);
        return 1.0f;
    }

    XrmDatabase db = XrmGetStringDatabase(res);
    XrmValue value;
    char *type = nullptr;

    if (XrmGetResource(db, "Xft.dpi", "Xft.Dpi", &type, &value))
    {
        float dpi = atof(value.addr);
        XCloseDisplay(dpy);
        return dpi / 96.0f;
    }

    XCloseDisplay(dpy);
    return 1.0f;

    // -----------------------------
    // LINUX (WAYLAND — fallback)
    // -----------------------------
#elif defined(__linux__)
    // Wayland does not give DPI easily; use scaling factor
    const char *wlScale = getenv("QT_SCALE_FACTOR");
    if (wlScale)
        return atof(wlScale);

    wlScale = getenv("GDK_SCALE"); // GNOME / GTK
    if (wlScale)
        return atof(wlScale);

    wlScale = getenv("XCURSOR_SIZE"); // fallback
    if (wlScale)
        return atof(wlScale) / 24.0f;

    return 1.0f;

#else
    return 1.0f;
#endif
}

void GraphicsRenderer::initialize(float width, float height)
{
    worldWidth = width;
    worldHeight = height;
    viewCenter = sf::Vector2f(worldWidth / 2.0f, worldHeight / 2.0f);
}

void GraphicsRenderer::handleEvent(sf::RenderWindow &window, const sf::Event &event, bool mouseCapture)
{
    // Mouse wheel zooming
    if (!mouseCapture)
    {
        if (const auto *wheel = event.getIf<sf::Event::MouseWheelScrolled>())
        {
            if (wheel->delta > 0)
                zoom = std::max(0.0001f, zoom * 0.9f); // zoom in (clamped)
            else
                zoom = std::min(10000.0f, zoom * 1.1f); // zoom out (clamped)
        }
    }

    // Window focus events
    if (event.is<sf::Event::FocusGained>())
    {
        isFocused = true;
    }
    else if (event.is<sf::Event::FocusLost>())
    {
        isFocused = false;
    }
}

void GraphicsRenderer::updatePanning(sf::RenderWindow &window, bool mouseCapture)
{
    // Mouse-based panning
    if (mouseCapture)
        return;

    sf::Vector2i mp = sf::Mouse::getPosition(window);
    sf::Vector2u size = window.getSize();

    bool insideWindow =
        mp.x >= 0 && mp.y >= 0 &&
        mp.x < (int)size.x && mp.y < (int)size.y;

    // When clicking down: decide if this drag is allowed
    if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
    {
        if (!isDragging)
        {
            // First frame of drag
            dragStartedInside = insideWindow;

            if (!dragStartedInside)
            {
                // Drag started from the title bar or outside window → ignore
                return;
            }

            lastMousePos = mp;
            isDragging = true;
        }

        if (!dragStartedInside)
            return; // ignore movement entirely

        // Normal panning
        sf::Vector2f worldLast = window.mapPixelToCoords(lastMousePos);
        sf::Vector2f worldNow = window.mapPixelToCoords(mp);

        viewCenter += (worldLast - worldNow);
        lastMousePos = mp;
    }
    else
    {
        // Reset when mouse released
        isDragging = false;
        dragStartedInside = false;
    }
}

void GraphicsRenderer::updateView(sf::RenderWindow &window)
{
    // Aspect-correct view
    sf::Vector2u windowSize = window.getSize();
    float windowAspect = static_cast<float>(windowSize.x) / windowSize.y;
    float worldAspect = worldWidth / worldHeight;

    sf::Vector2f viewSize;
    if (windowAspect >= worldAspect)
    {
        viewSize.y = worldHeight * zoom;
        viewSize.x = viewSize.y * windowAspect;
    }
    else
    {
        viewSize.x = worldWidth * zoom;
        viewSize.y = viewSize.x / windowAspect;
    }

    sf::View view;
    view.setCenter(viewCenter);
    view.setSize(sf::Vector2f(viewSize.x, -viewSize.y)); // negative y = y-up
    window.setView(view);
}

void GraphicsRenderer::drawGrid(sf::RenderWindow &window) const
{
    sf::View view = window.getView();
    sf::Vector2f viewSize = view.getSize();
    sf::Vector2f viewCenter = view.getCenter();

    // Compute visible world bounds and handle y-up (negative height)
    float left = viewCenter.x - viewSize.x / 2.f;
    float right = viewCenter.x + viewSize.x / 2.f;
    float top = viewCenter.y - viewSize.y / 2.f;
    float bottom = viewCenter.y + viewSize.y / 2.f;

    float xMin = std::min(left, right);
    float xMax = std::max(left, right);
    float yMin = std::min(top, bottom);
    float yMax = std::max(top, bottom);

    // Grid spacing in world units
    float gridSpacing = 1.0f; // adjust as needed

    // Compute line ranges
    float startX = std::floor(xMin / gridSpacing) * gridSpacing;
    float endX = std::ceil(xMax / gridSpacing) * gridSpacing;
    float startY = std::floor(yMin / gridSpacing) * gridSpacing;
    float endY = std::ceil(yMax / gridSpacing) * gridSpacing;

    sf::Color gridColor(130, 130, 130, 255);

    // Draw vertical lines
    for (float x = startX; x <= endX; x += gridSpacing)
    {
        sf::Vertex line[2];
        line[0].position = sf::Vector2f(x, yMin);
        line[0].color = gridColor;
        line[1].position = sf::Vector2f(x, yMax);
        line[1].color = gridColor;
        window.draw(line, 2, sf::PrimitiveType::Lines);
    }

    // Draw horizontal lines
    for (float y = startY; y <= endY; y += gridSpacing)
    {
        sf::Vertex line[2];
        line[0].position = sf::Vector2f(xMin, y);
        line[0].color = gridColor;
        line[1].position = sf::Vector2f(xMax, y);
        line[1].color = gridColor;
        window.draw(line, 2, sf::PrimitiveType::Lines);
    }
}

float GraphicsRenderer::getViewScale(const sf::RenderWindow &window) const
{
    // Get the current view to determine zoom level
    sf::View view = window.getView();
    sf::Vector2f viewSize = view.getSize();

    // Calculate scale based on view size relative to initial world size
    // The larger the view size, the more zoomed out we are
    float scale = std::abs(viewSize.x) / worldWidth;

    return scale;
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

    // Draw the main rectangle (quad)
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

    // Draw rounded endcaps
    const int circlePoints = 16; // number of segments for the circle
    sf::CircleShape cap(thickness * 0.5f, circlePoints);
    cap.setFillColor(color);

    // Center at start point
    cap.setPosition(a - sf::Vector2f(thickness * 0.5f, thickness * 0.5f));
    target.draw(cap);

    // Center at end point
    cap.setPosition(b - sf::Vector2f(thickness * 0.5f, thickness * 0.5f));
    target.draw(cap);
}

void GraphicsRenderer::drawSystem(sf::RenderWindow &window, const SpringSystem &system) const
{
    // Draw grid background first
    drawGrid(window);

    // Remember current world view to restore after drawing HUD (text)
    sf::View view = window.getView();

    // Get scale factor for zoom-independent sizing
    float viewScale = getViewScale(window);

    // Base sizes in pixels (what we want to see on screen)
    const float baseSpringThickness = 0.01f;
    const float baseNodeSize = 0.03f;
    const float baseArrowSize = 0.01f;

    // Convert to world units based on current zoom
    float springThickness = baseSpringThickness * viewScale;
    float nodeSize = baseNodeSize * viewScale;
    float arrowSize = baseArrowSize * viewScale;

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

        drawThickLine(window, p1, p2, springThickness, springColor);
    }

    // -------------------------
    // Draw nodes
    // -------------------------
    int index = 0;
    for (const auto &node : system.nodes)
    {
        sf::Vector2f pos(node.position[0] + system.displacement(index * 2),
                         node.position[1] + system.displacement(index * 2 + 1));

        if (node.constraint_type == Fixed)
        {
            sf::RectangleShape square(sf::Vector2f(nodeSize, nodeSize));
            square.setOrigin(sf::Vector2f(nodeSize / 2, nodeSize / 2));
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

            triangle.setPoint(0, sf::Vector2f(0, 0)); // tip
            triangle.setPoint(1, sf::Vector2f(-dx * nodeSize - dy * nodeSize / 2, -dy * nodeSize + dx * nodeSize / 2));
            triangle.setPoint(2, sf::Vector2f(-dx * nodeSize + dy * nodeSize / 2, -dy * nodeSize - dx * nodeSize / 2));

            triangle.setPosition(pos);
            triangle.setFillColor(sf::Color::Yellow);
            window.draw(triangle);
        }
        else
        {
            float radius = nodeSize / 2.0f;
            sf::CircleShape circle(radius);
            circle.setOrigin(sf::Vector2f(radius, radius));
            circle.setPosition(pos);
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);
        }

        // Draw text label in world coordinates
        sf::Text labelText(font);
        labelText.setString(std::to_string(index + 1));
        labelText.setCharacterSize(25); // Scale with zoom
        labelText.setStyle(sf::Text::Bold);
        labelText.setFillColor(sf::Color::White);
        labelText.setOutlineColor(sf::Color::Black);
        labelText.setOutlineThickness(4.f);

        // Center the text
        sf::FloatRect lb = labelText.getLocalBounds();
        sf::Vector2f center(lb.position.x + lb.size.x / 2.f, lb.position.y + lb.size.y / 2.f);
        labelText.setOrigin(center);

        // Position in world coordinates
        labelText.setPosition(pos);
        labelText.setScale(sf::Vector2f(viewScale / 1000, -viewScale / 1000)); // Flip text to match y-up coordinate system

        window.draw(labelText);

        ++index;
    }
    // -------------------------
    // Draw forces as arrows
    // -------------------------
    for (int i = 0; i < system.nodes.size(); ++i)
    {
        const Node &node = system.nodes[i];

        // scale factor for force visualization (increase to make arrows shorter)
        const float arrowScale = 500.0f;
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

            sf::Vertex arrow[3];
            arrow[0].position = end;
            arrow[0].color = sf::Color::Magenta;
            arrow[1].position = end - unit * arrowSize + perp * (arrowSize / 2);
            arrow[1].color = sf::Color::Magenta;
            arrow[2].position = end - unit * arrowSize - perp * (arrowSize / 2);
            arrow[2].color = sf::Color::Magenta;

            window.draw(arrow, 3, sf::PrimitiveType::Triangles);
        }
    }
}
