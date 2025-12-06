#include "graphics.h"

GraphicsRenderer::GraphicsRenderer(FEMSystem const &springSystem)
    : worldWidth(1.2f),
      worldHeight(1.2f),
      viewCenter(worldWidth / 2.0f, worldHeight / 2.0f),
      zoom(20.0f),
      isDragging(false),
      isFocused(true),
      dragStartedInside(false),
      font(),
      text(font),
      system(springSystem)
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
    // WINDOWS
    // -----------------------------
#if defined(_WIN32)
    UINT dpi = GetDpiForWindow(handle);
    return dpi / 96.0f;

    // -----------------------------
    // MACOS
    // -----------------------------
#elif defined(__APPLE__)
    // macOS: Get the window's screen and its backing scale factor
    @autoreleasepool
    {
        NSWindow *nsWindow = (__bridge NSWindow *)handle;
        if (nsWindow && nsWindow.screen)
        {
            CGFloat scale = nsWindow.screen.backingScaleFactor;
            return static_cast<float>(scale);
        }
        // Fallback to main screen
        CGFloat scale = [[NSScreen mainScreen] backingScaleFactor];
        return static_cast<float>(scale);
    }

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
    if (!db)
    {
        XCloseDisplay(dpy);
        return 1.0f;
    }

    XrmValue value;
    char *type = nullptr;
    float dpiScale = 1.0f;

    if (XrmGetResource(db, "Xft.dpi", "Xft.Dpi", &type, &value))
    {
        if (value.addr)
        {
            float dpi = atof(value.addr);
            if (dpi > 0.0f)
            {
                dpiScale = dpi / 96.0f;
            }
        }
    }

    XrmDestroyDatabase(db);
    XCloseDisplay(dpy);
    return dpiScale;

    // -----------------------------
    // LINUX (WAYLAND — fallback)
    // -----------------------------
#elif defined(__linux__)
    // Wayland: Check various scaling environment variables
    const char *wlScale = nullptr;

    // Try GDK_SCALE first (GNOME/GTK - integer scale factor)
    wlScale = getenv("GDK_SCALE");
    if (wlScale && *wlScale)
    {
        float scale = atof(wlScale);
        if (scale > 0.0f)
            return scale;
    }

    // Try QT_SCALE_FACTOR (Qt applications)
    wlScale = getenv("QT_SCALE_FACTOR");
    if (wlScale && *wlScale)
    {
        float scale = atof(wlScale);
        if (scale > 0.0f)
            return scale;
    }

    // Try GDK_DPI_SCALE (fractional scaling)
    wlScale = getenv("GDK_DPI_SCALE");
    if (wlScale && *wlScale)
    {
        float scale = atof(wlScale);
        if (scale > 0.0f)
            return scale;
    }

    // XCURSOR_SIZE is not reliable for DPI scaling
    // Default cursor size is 24, but this varies by theme

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

void GraphicsRenderer::drawCubicBezierThick(sf::RenderTarget &target,
                                            const sf::Vector2f &p0,
                                            const sf::Vector2f &p1,
                                            const sf::Vector2f &p2,
                                            const sf::Vector2f &p3,
                                            float thickness,
                                            const sf::Color &color,
                                            int segments) const
{
    if (segments < 1)
        segments = 1;

    sf::Vector2f previousPoint = p0;

    for (int i = 1; i <= segments; ++i)
    {
        float t = static_cast<float>(i) / segments;
        float u = 1.0f - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;

        // Standard Cubic Bezier Formula:
        // B(t) = (1-t)^3*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3
        sf::Vector2f currentPoint =
            (uuu * p0) +
            (3 * uu * t * p1) +
            (3 * u * tt * p2) +
            (ttt * p3);

        // Reuse existing thick line drawer for the segment
        drawThickLine(target, previousPoint, currentPoint, thickness, color);
        previousPoint = currentPoint;
    }
}

void GraphicsRenderer::drawSystem(sf::RenderWindow &window) const
{
    // Draw grid background first
    drawGrid(window);

    // Remember current world view to restore after drawing HUD (text)
    sf::View view = window.getView();

    // Get scale factor for zoom-independent sizing
    float viewScale = getViewScale(window);

    // Base sizes in pixels (what we want to see on screen)
    const float baseBeamThickness = 0.01f;
    const float baseNodeSize = 0.03f;
    const float baseArrowSize = 0.01f;
    const int curveSegments = 24;

    // Convert to world units based on current zoom
    float beamThickness = baseBeamThickness * viewScale;
    float nodeSize = baseNodeSize * viewScale;
    float arrowSize = baseArrowSize * viewScale;

    // -------------------------
    // Draw beams with stress-based colors
    // -------------------------
    for (size_t i = 0; i < system.beams.size(); ++i)
    {
        const auto &beam = system.beams[i];
        int n1_idx = beam.nodes[0];
        int n2_idx = beam.nodes[1];

        sf::Color beamColor = getStressColor(beam.stress, system.min_stress, system.max_stress);

        // 1. Get Displaced Endpoints positions (P0 and P3 for Bezier)
        sf::Vector2f p0_displaced(
            system.nodes[n1_idx].position[0] + system.displacement(n1_idx * 3),
            system.nodes[n1_idx].position[1] + system.displacement(n1_idx * 3 + 1));

        sf::Vector2f p3_displaced(
            system.nodes[n2_idx].position[0] + system.displacement(n2_idx * 3),
            system.nodes[n2_idx].position[1] + system.displacement(n2_idx * 3 + 1));

        // NEW: --- Calculate Bezier Control Points based on Rotation ---

        // Get the rotational displacement (theta) at each node.
        // IMPORTANT: This assumes your system is solving 3 DOF per node (x, y, theta).
        // Theta is usually at index 3*i + 2.
        float theta1_rad = static_cast<float>(system.displacement(n1_idx * 3 + 2));
        float theta2_rad = static_cast<float>(system.displacement(n2_idx * 3 + 2));

        // Calculate the vector of the displaced beam chord
        sf::Vector2f chordDir = p3_displaced - p0_displaced;
        float chordLength = std::sqrt(chordDir.x * chordDir.x + chordDir.y * chordDir.y);

        // Avoid division by zero for zero-length beams
        if (chordLength < 1e-6f)
        {
            drawThickLine(window, p0_displaced, p3_displaced, beamThickness, beamColor);
            continue;
        }

        // Current angle of the beam's chord in global space
        float currentChordAngle = std::atan2(chordDir.y, chordDir.x);

        // The actual tangent angle at the node is the chord angle PLUS the nodal rotation.
        // Note: In some formulations, theta is absolute. If your visualization looks wrong,
        // try removing currentChordAngle and just using theta1_rad.
        // Based on typical frame formulations, theta is relative to the beam axis. Let's try relative first.
        // Wait, no, in standard global stiffness method, theta is the absolute global rotation of the node. Let's use absolute.

        // REVISION: In standard 2D frame analysis, the calculated Theta value in the solution vector
        // is usually the *change* in rotation from the initial state (0).
        // We need the initial angle of the undeformed beam.
        sf::Vector2f initialDir(
            static_cast<float>(system.nodes[n2_idx].position[0] - system.nodes[n1_idx].position[0]),
            static_cast<float>(system.nodes[n2_idx].position[1] - system.nodes[n1_idx].position[1]));
        float initialAngle = std::atan2(initialDir.y, initialDir.x);

        // Final tangent angles = Initial Angle + Calculated Rotation Change
        float tangent1Angle = initialAngle + theta1_rad;
        float tangent2Angle = initialAngle + theta2_rad;

        // Calculate control point distance. A heuristic of L/3 or L/4 works well visually.
        float controlDist = chordLength * 0.33f;

        // Control Point 1 (P1): Start at P0, move along tangent 1 direction
        sf::Vector2f p1_control = p0_displaced + sf::Vector2f(
                                                     std::cos(tangent1Angle) * controlDist,
                                                     std::sin(tangent1Angle) * controlDist);

        // Control Point 2 (P2): Start at P3, move BACKWARDS along tangent 2 direction
        sf::Vector2f p2_control = p3_displaced - sf::Vector2f(
                                                     std::cos(tangent2Angle) * controlDist,
                                                     std::sin(tangent2Angle) * controlDist);

        // Draw the curved beam
        drawCubicBezierThick(window, p0_displaced, p1_control, p2_control, p3_displaced, beamThickness, beamColor, curveSegments);

        // ---------------------------------------------------------

        // Draw beam number label at the center of the chord (could be improved to follow curve, but this is okay)
        sf::Vector2f center = (p0_displaced + p3_displaced) * 0.5f;
        sf::Text springLabel(font);
        springLabel.setString(std::to_string(i));
        springLabel.setCharacterSize(25);
        springLabel.setStyle(sf::Text::Regular);
        springLabel.setFillColor(sf::Color::Black);
        springLabel.setOutlineColor(sf::Color::White);
        springLabel.setOutlineThickness(4.f * viewScale / 20.0f); // Scale outline too

        sf::FloatRect lb = springLabel.getLocalBounds();
        sf::Vector2f labelCenter(lb.getCenter());
        springLabel.setOrigin(labelCenter);

        springLabel.setPosition(center);
        // Adjust scale to keep text readable relative to zoom
        float textScale = viewScale / 700.0f;
        springLabel.setScale(sf::Vector2f(textScale, -textScale));

        window.draw(springLabel);
    }

    // -------------------------
    // Draw nodes
    // -------------------------
    int index = 0;
    for (const auto &node : system.nodes)
    {
        sf::Vector2f pos(node.position[0] + system.displacement(index * 3),
                         node.position[1] + system.displacement(index * 3 + 1));

        if (node.constraint_type == FixedPin)
        {
            // Draw a thick red "X" for FixedPin using drawThickLine
            float halfSize = nodeSize / 2.0f;
            float thickness = nodeSize * 0.4f; // Adjust thickness as needed

            sf::Vector2f p1 = pos + sf::Vector2f(-halfSize, -halfSize);
            sf::Vector2f p2 = pos + sf::Vector2f(halfSize, halfSize);
            sf::Vector2f p3 = pos + sf::Vector2f(-halfSize, halfSize);
            sf::Vector2f p4 = pos + sf::Vector2f(halfSize, -halfSize);

            drawThickLine(window, p1, p2, thickness, sf::Color::Red);
            drawThickLine(window, p3, p4, thickness, sf::Color::Red);
        }
        else if (node.constraint_type == Fixed)
        {
            float halfSize = nodeSize / 2.0f;
            sf::RectangleShape square(sf::Vector2f(nodeSize, nodeSize));
            square.setOrigin(sf::Vector2f(halfSize, halfSize));
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
        float fx = system.forces(i * 3) / arrowScale;
        float fy = system.forces(i * 3 + 1) / arrowScale;

        if (std::abs(fx) < 1e-6f && std::abs(fy) < 1e-6f)
            continue;

        sf::Vector2f start(node.position[0] + system.displacement(i * 3),
                           node.position[1] + system.displacement(i * 3 + 1));
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

void GraphicsRenderer::centerView()
{
    if (system.nodes.empty())
        return;

    sf::Vector2f sum(0.f, 0.f);
    for (const auto &node : system.nodes)
    {
        sum.x += node.position[0];
        sum.y += node.position[1];
    }
    viewCenter = sf::Vector2f(sum.x / system.nodes.size(), sum.y / system.nodes.size());
}

void GraphicsRenderer::autoZoomToFit()
{
    if (system.nodes.empty())
        return;

    float minX = system.nodes[0].position[0];
    float maxX = system.nodes[0].position[0];
    float minY = system.nodes[0].position[1];
    float maxY = system.nodes[0].position[1];

    for (const auto &node : system.nodes)
    {
        if (node.position[0] < minX)
            minX = node.position[0];
        if (node.position[0] > maxX)
            maxX = node.position[0];
        if (node.position[1] < minY)
            minY = node.position[1];
        if (node.position[1] > maxY)
            maxY = node.position[1];
    }

    float padding = 0.1f; // Add some padding around the system
    float width = (maxX - minX) + padding * 2;
    float height = (maxY - minY) + padding * 2;

    // Center view
    viewCenter = sf::Vector2f((minX + maxX) / 2.0f, (minY + maxY) / 2.0f);

    // Adjust zoom to fit
    float aspectRatio = worldWidth / worldHeight;
    if (width / height > aspectRatio)
    {
        zoom = width / worldWidth;
    }
    else
    {
        zoom = height / worldHeight;
    }
}