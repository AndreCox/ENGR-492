#include "gui_handler.h"
#include <imgui-SFML.h>
#include <fstream>
#include <cstdint>

GUIHandler::GUIHandler(SpringSystem &system, GraphicsRenderer &renderer, sf::RenderWindow &window)
    : window(window), spring_system(system), renderer(renderer), show_system_controls(true), show_node_editor(false), show_spring_editor(false)
{
    ImGuiIO &io = ImGui::GetIO();
    float dpiScale = renderer.GetDPIScale(window.getNativeHandle());
    io.FontGlobalScale = dpiScale;             // Scale text
    ImGui::GetStyle().ScaleAllSizes(dpiScale); // Scale widgets
}

void GUIHandler::processEvent(const sf::Event &event)
{
    ImGui::SFML::ProcessEvent(window, event);
}

void GUIHandler::render()
{
    ImGui::SFML::Update(window, deltaClock.restart());

    systemControls();
    nodeEditor();
    springEditor();
    handleSavePopup();
    handleLoadPopup();
    headerBar();
}

void GUIHandler::nodeEditor()
{
    if (!show_node_editor)
        return;

    ImGui::Begin("Node Editor", &show_node_editor, ImGuiWindowFlags_AlwaysAutoResize);

    ImGui::Text("Nodes in System:");
    ImGui::Separator();

    bool constraints_changed = false;

    for (int i = 0; i < spring_system.nodes.size(); ++i)
    {
        ImGui::PushID(i);

        const auto &node = spring_system.nodes[i];

        ImGui::Text("Node %d", i + 1);

        // Editable reference to the node
        auto &editable_node = spring_system.nodes[i];

        // Ensure the combo items match the enum mapping used below (Free=0, Fixed=1, FixedPin=2, Slider=3)
        const char *constraint_items[] = {"Free", "Fixed", "Fixed Pin", "Slider"};
        int current_constraint = 0;
        switch (editable_node.constraint_type)
        {
        case Free:
            current_constraint = 0;
            break;
        case Fixed:
            current_constraint = 1;
            break;
        case FixedPin:
            current_constraint = 2;
            break;
        case Slider:
            current_constraint = 3;
            break;
        default:
            current_constraint = 0;
            break;
        }

        if (ImGui::Combo("Constraint", &current_constraint, constraint_items, IM_ARRAYSIZE(constraint_items)))
        {
            switch (current_constraint)
            {
            case 0:
                editable_node.constraint_type = Free;
                break;
            case 1:
                editable_node.constraint_type = Fixed;
                break;
            case 2:
                editable_node.constraint_type = FixedPin;
                break;
            case 3:
                editable_node.constraint_type = Slider;
                break;
            }
            constraints_changed = true;
        }

        // If Slider, allow angle adjustment
        if (editable_node.constraint_type == Slider)
        {
            // Show and update slider angle if node is a Slider
            float angle_deg = editable_node.constraint_angle; // Node must be constructed as: Node(x, y, Slider, angle_deg)
            if (ImGui::SliderFloat("Slider Angle (deg)", &angle_deg, 0.0f, 360.0f))
            {
                editable_node.constraint_angle = angle_deg;
                constraints_changed = true;
            }
        }

        // Remove node button
        ImGui::SameLine();
        if (ImGui::Button("Remove Node"))
        {
            const int removed_index = i;

            // Remove any springs that reference this node
            for (int s = 0; s < static_cast<int>(spring_system.springs.size()); ++s)
            {
                auto &sp = spring_system.springs[s];
                if (sp.nodes[0] == removed_index || sp.nodes[1] == removed_index)
                {
                    spring_system.springs.erase(spring_system.springs.begin() + s);
                    --s; // step back after erase
                }
            }

            // Decrement node indices in remaining springs that were after the removed node
            for (auto &sp : spring_system.springs)
            {
                if (sp.nodes[0] > removed_index)
                    --sp.nodes[0];
                if (sp.nodes[1] > removed_index)
                    --sp.nodes[1];
            }

            // Erase the node
            spring_system.nodes.erase(spring_system.nodes.begin() + removed_index);

            ImGui::PopID();

            // Rebuild/solve system after structural change
            spring_system.solve_system();

            // Adjust loop index to account for erased element
            --i;
            continue;
        }

        ImGui::Separator();
        ImGui::PopID();
    }

    // --- New node creation UI ---
    ImGui::Separator();
    ImGui::Text("Create New Node:");
    static float new_x = 0.0f;
    static float new_y = 0.0f;
    static int new_constraint = 0;
    static float new_angle = 0.0f;
    const char *constraint_items[] = {"Free", "Fixed", "Fixed Pin", "Slider"};

    ImGui::InputFloat("X (m)", &new_x);
    ImGui::InputFloat("Y (m)", &new_y);
    ImGui::Combo("Constraint##NewNode", &new_constraint, constraint_items, IM_ARRAYSIZE(constraint_items));
    if (new_constraint == 3) // Slider
    {
        ImGui::SliderFloat("Slider Angle (deg)##NewNode", &new_angle, 0.0f, 360.0f);
    }

    ImGui::SameLine();
    if (ImGui::Button("Add Node"))
    {
        // Map int to enum and construct node. Adjust constructor call if your Node type differs.
        switch (new_constraint)
        {
        case 0:
            spring_system.nodes.emplace_back(new_x, new_y, Free, 0.0f);
            break;
        case 1:
            spring_system.nodes.emplace_back(new_x, new_y, Fixed, 0.0f);
            break;
        case 2:
            spring_system.nodes.emplace_back(new_x, new_y, FixedPin, 0.0f);
            break;
        case 3:
            spring_system.nodes.emplace_back(new_x, new_y, Slider, new_angle);
            break;
        default:
            spring_system.nodes.emplace_back(new_x, new_y, Free, 0.0f);
            break;
        }

        // Optionally bump the X a bit so successive adds are easier
        new_x += 0.5f;

        // Re-solve so system matrices/vectors are rebuilt (solve_system should handle resizing state)
        spring_system.solve_system();
    }

    ImGui::End();

    if (constraints_changed)
    {
        spring_system.solve_system();
    }
}

void GUIHandler::springEditor()
{
    if (!show_spring_editor)
        return;

    ImGui::Begin("Spring Editor", &show_spring_editor, ImGuiWindowFlags_AlwaysAutoResize);

    // Prepare node labels for combos
    std::vector<std::string> node_labels;
    node_labels.reserve(spring_system.nodes.size());
    for (size_t i = 0; i < spring_system.nodes.size(); ++i)
        node_labels.push_back("Node " + std::to_string(i + 1));
    std::vector<const char *> node_items;
    node_items.reserve(node_labels.size());
    for (auto &s : node_labels)
        node_items.push_back(s.c_str());

    bool springs_changed = false;

    // List existing springs and allow changing endpoints or removing
    for (int i = 0; i < static_cast<int>(spring_system.springs.size()); ++i)
    {
        ImGui::PushID(i);
        auto &spring = spring_system.springs[i];

        ImGui::Text("Spring %d", i + 1);

        // Node A combo
        int node_a = spring.nodes[0]; // adjust field name if your Spring uses a different member
        if (ImGui::Combo("Node A", &node_a, node_items.data(), static_cast<int>(node_items.size())))
        {
            spring.nodes[0] = node_a;
            springs_changed = true;
        }

        // Node B combo
        int node_b = spring.nodes[1]; // adjust field name if your Spring uses a different member
        if (ImGui::Combo("Node B", &node_b, node_items.data(), static_cast<int>(node_items.size())))
        {
            spring.nodes[1] = node_b;
            springs_changed = true;
        }

        // Optional: display stress if available
        ImGui::SameLine();
        ImGui::Text("Stress: %.2f", spring.stress);

        // Remove button
        if (ImGui::Button("Remove Spring"))
        {
            spring_system.springs.erase(spring_system.springs.begin() + i);
            springs_changed = true;
            ImGui::PopID();
            --i; // step back to account for erased element
            continue;
        }

        ImGui::Separator();
        ImGui::PopID();
    }

    // Add new spring controls
    static int new_node_a = 0;
    static int new_node_b = 0;
    if (!node_items.empty())
    {
        ImGui::Text("Create New Spring:");
        ImGui::Combo("New Node A", &new_node_a, node_items.data(), static_cast<int>(node_items.size()));
        ImGui::Combo("New Node B", &new_node_b, node_items.data(), static_cast<int>(node_items.size()));

        ImGui::SameLine();
        if (ImGui::Button("Add Spring"))
        {
            if (new_node_a != new_node_b)
            {
                // Try to construct a new spring using a common constructor; adjust if your Spring type differs
                // default area (m^2) and Young's modulus (Pa) — adjust defaults or expose as UI controls elsewhere
                static float new_area = 0.1963; // cross-sectional area
                static float new_E = 30e6f;     // Young's modulus (e.g. 200 GPa)
                BeamProperties properties;
                properties.area = new_area;
                properties.material.youngs_modulus = new_E;
                // if (new_area <= 0.0f)
                //     new_area = 1.0f;
                // if (new_E <= 0.0f)
                //     new_E = 1.0f;

                // construct spring with node indices a/b and material properties A and E
                spring_system.springs.emplace_back(new_node_a, new_node_b, properties);
                springs_changed = true;
            }
        }
    }
    else
    {
        ImGui::TextDisabled("No nodes available to create springs.");
    }

    ImGui::End();

    if (springs_changed)
    {
        spring_system.solve_system();
    }
}

void GUIHandler::systemControls()
{

    if (!show_system_controls)
        return;

    ImGui::Begin("System Controls", &show_system_controls, ImGuiWindowFlags_AlwaysAutoResize);

    ImGui::Text("2D Truss System Solver");
    ImGui::Separator();

    bool forces_changed = false;
    for (int i = 0; i < spring_system.nodes.size(); ++i)
    {
        if (spring_system.nodes[i].constraint_type == Free || spring_system.nodes[i].constraint_type == Slider)
        {
            float fx = spring_system.forces(i * 2);
            float fy = spring_system.forces(i * 2 + 1);
            ImGui::PushID(i);

            ImGui::Text("Node %d Forces:", i + 1);

            // Fx controls
            ImGui::Text("Fx:");
            ImGui::SameLine();
            if (ImGui::SliderFloat("##FxSlider", &fx, -10000.0f, 10000.0f))
            {
                spring_system.forces(i * 2) = fx;
                forces_changed = true;
            }
            ImGui::SameLine();
            if (ImGui::InputFloat("##FxInput", &fx))
            {
                spring_system.forces(i * 2) = fx;
                forces_changed = true;
            }

            // Fy controls
            ImGui::Text("Fy:");
            ImGui::SameLine();
            if (ImGui::SliderFloat("##FySlider", &fy, -10000.0f, 10000.0f))
            {
                spring_system.forces(i * 2 + 1) = fy;
                forces_changed = true;
            }
            ImGui::SameLine();
            if (ImGui::InputFloat("##FyInput", &fy))
            {
                spring_system.forces(i * 2 + 1) = fy;
                forces_changed = true;
            }

            ImGui::Separator();
            ImGui::PopID();
        }
    }

    if (forces_changed)
        spring_system.solve_system();

    // Solution display
    ImGui::Text("Solution:");
    for (int i = 0; i < spring_system.nodes.size(); ++i)
    {
        ImGui::Text("Node %d: u=%.6f m, v=%.6f m", i + 1,
                    spring_system.displacement(i * 2),
                    spring_system.displacement(i * 2 + 1));
    }

    // Spring stresses
    ImGui::Text("Spring Stresses (MPa):");
    int index = 1;
    for (const auto &spring : spring_system.springs)
    {
        sf::Color color = renderer.getStressColor(spring.stress, spring_system.min_stress, spring_system.max_stress);
        ImGui::TextColored(ImVec4(color.r / 255.f, color.g / 255.f, color.b / 255.f, 1.f),
                           "Spring %d: %.2f MPa", index++, spring.stress);
    }

    ImGui::End();
}

void GUIHandler::handleLoadPopup()
{
    static char filename_buf[512] = "system.bin";
    static bool trigger_load_read = false;

    if (request_load_popup)
    {
        ImGui::OpenPopup("Load From File");
        request_load_popup = false;
    }

    if (ImGui::BeginPopupModal("Load From File", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::InputText("Filename", filename_buf, sizeof(filename_buf));

        if (ImGui::Button("Load"))
        {
            trigger_load_read = true;
            ImGui::CloseCurrentPopup();
        }

        ImGui::SameLine();
        if (ImGui::Button("Cancel"))
        {
            trigger_load_read = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    if (!trigger_load_read)
        return;

    // consume the trigger
    trigger_load_read = false;

    std::ifstream ifs(filename_buf, std::ios::binary);
    if (!ifs)
    {
        ImGui::OpenPopup("Load Error");
        return;
    }

    // Try to read nodes
    try
    {
        uint32_t node_count = 0;
        ifs.read(reinterpret_cast<char *>(&node_count), sizeof(node_count));
        if (!ifs)
            throw std::runtime_error("read node count");

        // Temporary containers
        std::vector<std::tuple<float, float, int32_t, float>> nodes_data;
        nodes_data.reserve(node_count);

        for (uint32_t i = 0; i < node_count; ++i)
        {
            float x = 0.f, y = 0.f;
            int32_t type = 0;
            float angle = 0.f;

            ifs.read(reinterpret_cast<char *>(&x), sizeof(x));
            ifs.read(reinterpret_cast<char *>(&y), sizeof(y));
            ifs.read(reinterpret_cast<char *>(&type), sizeof(type));
            ifs.read(reinterpret_cast<char *>(&angle), sizeof(angle));

            if (!ifs)
                throw std::runtime_error("read node entry");

            nodes_data.emplace_back(x, y, type, angle);
        }

        uint32_t spring_count = 0;
        ifs.read(reinterpret_cast<char *>(&spring_count), sizeof(spring_count));
        if (!ifs)
            throw std::runtime_error("read spring count");

        std::vector<std::tuple<int32_t, int32_t, float>> springs_data;
        springs_data.reserve(spring_count);

        for (uint32_t i = 0; i < spring_count; ++i)
        {
            int32_t a = 0, b = 0;
            float stress = 0.f;
            ifs.read(reinterpret_cast<char *>(&a), sizeof(a));
            ifs.read(reinterpret_cast<char *>(&b), sizeof(b));
            ifs.read(reinterpret_cast<char *>(&stress), sizeof(stress));

            if (!ifs)
                throw std::runtime_error("read spring entry");

            springs_data.emplace_back(a, b, stress);
        }

        uint32_t fcount = 0;
        ifs.read(reinterpret_cast<char *>(&fcount), sizeof(fcount));
        if (!ifs)
            throw std::runtime_error("read forces count");

        std::vector<float> forces_data;
        forces_data.resize(fcount ? fcount : 0);
        if (fcount > 0)
        {
            ifs.read(reinterpret_cast<char *>(forces_data.data()), sizeof(float) * fcount);
            if (!ifs)
                throw std::runtime_error("read forces");
        }

        // If we reached here, data read successfully. Replace system content.
        spring_system.nodes.clear();
        spring_system.nodes.reserve(nodes_data.size());
        for (auto &t : nodes_data)
        {
            float x, y, angle;
            int32_t type;
            std::tie(x, y, type, angle) = t;

            // Construct node directly with saved parameters (Node(float x, float y, ConstraintType ct, float angle))
            spring_system.nodes.emplace_back(x, y, static_cast<ConstraintType>(type), angle);
        }

        spring_system.springs.clear();
        spring_system.springs.reserve(springs_data.size());
        for (auto &t : springs_data)
        {
            int32_t a, b;
            float stress;
            std::tie(a, b, stress) = t;

            // Construct spring with the saved node indices and reasonable default A and E,
            // then restore the saved stress value.
            // TODO: Fix spring loading
            // spring_system.springs.emplace_back(a, b, );
            auto &s = spring_system.springs.back();
            s.stress = stress;
        }

        // Restore forces
        if (forces_data.size() == spring_system.forces.size())
        {
            for (size_t i = 0; i < forces_data.size(); ++i)
            {
                spring_system.forces(i) = forces_data[i];
            }
        }

        // Recompute system structures/state
        spring_system.solve_system();
    }
    catch (...)
    {
        ImGui::OpenPopup("Load Error");
    }
}

void GUIHandler::handleSavePopup()
{
    static char filename_buf[512] = "system.bin";
    static bool trigger_save_write = false;

    if (request_save_popup) // Check the flag
    {
        ImGui::OpenPopup("Save As"); // Open it here
        request_save_popup = false;
    }

    if (ImGui::BeginPopupModal("Save As", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::InputText("Filename", filename_buf, sizeof(filename_buf));

        if (ImGui::Button("Save"))
        {
            trigger_save_write = true;
            ImGui::CloseCurrentPopup();
        }

        ImGui::SameLine();
        if (ImGui::Button("Cancel"))
        {
            trigger_save_write = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    if (trigger_save_write)
    {
        trigger_save_write = false;

        std::ofstream ofs(filename_buf, std::ios::binary);

        if (!ofs)
        {
            ImGui::OpenPopup("Save Error");
        }
        else
        {
            // Write nodes
            uint32_t node_count = spring_system.nodes.size();
            ofs.write((char *)&node_count, sizeof(node_count));
            for (auto &n : spring_system.nodes)
            {
                float x = n.position[0];
                float y = n.position[1];
                ofs.write((char *)&x, sizeof(float));
                ofs.write((char *)&y, sizeof(float));

                int32_t type = (int32_t)n.constraint_type;
                ofs.write((char *)&type, sizeof(type));

                float angle = n.constraint_angle;
                ofs.write((char *)&angle, sizeof(angle));
            }

            // Write springs
            uint32_t spring_count = spring_system.springs.size();
            ofs.write((char *)&spring_count, sizeof(spring_count));
            for (auto &s : spring_system.springs)
            {
                int32_t a = s.nodes[0];
                int32_t b = s.nodes[1];
                ofs.write((char *)&a, sizeof(a));
                ofs.write((char *)&b, sizeof(b));

                float stress = s.stress;
                ofs.write((char *)&stress, sizeof(stress));
            }

            // Write forces
            uint32_t fcount = spring_system.forces.size();
            ofs.write((char *)&fcount, sizeof(fcount));
            if (fcount > 0)
                ofs.write((char *)spring_system.forces.data(),
                          sizeof(float) * fcount);
        }
    }
}

void GUIHandler::materialEditor()
{
    if (!show_material_editor)
        return;

    ImGui::Begin("Material Editor", &show_material_editor, ImGuiWindowFlags_AlwaysAutoResize);

    // Material editor content goes here

    ImGui::End();
}

void GUIHandler::headerBar()
{
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("New System", "Ctrl+N"))
            {
                spring_system.nodes.clear();
                spring_system.springs.clear();
                spring_system.solve_system();
            }
            if (ImGui::MenuItem("Open...", "Ctrl+O"))
            {
                request_load_popup = true;
            }

            if (ImGui::MenuItem("Save", "Ctrl+S"))
            {
                request_save_popup = true;
            }

            ImGui::Separator();
            if (ImGui::MenuItem("Exit", "Alt+F4"))
            {
                window.close();
            }
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("View"))
        {
            if (ImGui::MenuItem("Auto View"))
            {
                renderer.centerView();
                renderer.autoZoomToFit();
            }
            if (ImGui::MenuItem("Center View"))
            {
                renderer.centerView();
            }
            if (ImGui::MenuItem("Auto Zoom"))
            {
                renderer.autoZoomToFit();
            }
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("System"))
        {
            if (ImGui::MenuItem("System Controls"))
            {
                show_system_controls = !show_system_controls;
            }
            if (ImGui::MenuItem("System Node Editor"))
            {
                show_node_editor = !show_node_editor;
            }
            if (ImGui::MenuItem("System Spring Editor"))
            {
                show_spring_editor = !show_spring_editor;
            }

            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("Materials & Properties"))
        {
            if (ImGui::MenuItem("Material Editor"))
            {
                show_material_editor = !show_material_editor;
            }
            if (ImGui::MenuItem("Profile Editor"))
            {
                show_profile_editor = !show_profile_editor;
            }
            if (ImGui::MenuItem("Section Editor"))
            {
                show_section_editor = !show_section_editor;
            }

            ImGui::EndMenu();
        }

        float frame_rate = ImGui::GetIO().Framerate;
        ImGui::SameLine(ImGui::GetWindowWidth() - 100);
        ImGui::Text("FPS: %.1f", frame_rate);

        ImGui::EndMainMenuBar();
    }
}