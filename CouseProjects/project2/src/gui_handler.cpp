#include "gui_handler.h"
#include <imgui-SFML.h>
#include <fstream>
#include <cstdint>
#include "serialization.h"

GUIHandler::GUIHandler(FEMSystem &system, GraphicsRenderer &renderer, sf::RenderWindow &window)
    : window(window), fem_system(system), renderer(renderer),
      show_system_controls(true), show_node_editor(false), show_beam_editor(false)
{

    float detected = 1.0f;
    try
    {
        detected = renderer.GetDPIScale(window.getNativeHandle());
    }
    catch (...)
    {
        // fallback to 1.0
        detected = 1.0f;
    }

    applyDPIScale(detected);
}

void GUIHandler::reloadFontsForDPI(float dpiScale)
{
    ImGuiIO &io = ImGui::GetIO();
    float font_size = 16.0f * current_dpi_scale;
    io.Fonts->Clear();

    ImFont *roboto = io.Fonts->AddFontFromFileTTF("resources/fonts/Roboto-Regular.ttf", font_size);
    if (!roboto)
    {
        std::cerr << "Failed to load Roboto\n";
        // fallback to default font
        io.Fonts->AddFontDefault();
    }
    io.Fonts->Build();
    (void)ImGui::SFML::UpdateFontTexture();
}

void GUIHandler::applyDPIScale(float scale)
{
    ImGuiIO &io = ImGui::GetIO();

    // Reset style and font to base, then apply new scale exactly once.
    ImGui::GetStyle() = base_imgui_style;
    ImGui::GetStyle().ScaleAllSizes(scale);

    current_dpi_scale = scale;

    reloadFontsForDPI(scale);

    io.FontGlobalScale = base_font_global_scale;

    std::cout << "DPI Scale set to: " << scale << " (" << (scale * 100.0f) << "%)\n";
}

void GUIHandler::processEvent(const sf::Event &event)
{
    ImGui::SFML::ProcessEvent(window, event);

    if (ImGui::GetIO().WantCaptureKeyboard)
        return;

    using sc = sf::Keyboard::Scancode;

    // Key pressed event
    if (event.is<sf::Event::KeyPressed>())
    {
        const auto &e = event.getIf<sf::Event::KeyPressed>();

        bool ctrl = sf::Keyboard::isKeyPressed(sc::LControl) ||
                    sf::Keyboard::isKeyPressed(sc::RControl);

        sf::Keyboard::Scancode code = e->scancode;

        if (ctrl && code == sc::S)
        {
            request_save_popup = true;
        }
        else if (ctrl && code == sc::O)
        {
            request_load_popup = true;
        }
        else if (ctrl && code == sc::N)
        {
            fem_system.nodes.clear();
            fem_system.beams.clear();
            fem_system.solve_system();
        }
    }
}

void GUIHandler::render()
{
    if (pending_dpi_scale > 0.0f)
    {
        applyDPIScale(pending_dpi_scale);
        pending_dpi_scale = 0.0f; // Reset the flag
    }

    ImGui::SFML::Update(window, deltaClock.restart());

    systemControls();
    nodeEditor();
    beamEditor();
    materialEditor();
    profileEditor();
    handleSavePopup();
    handleLoadPopup();
    handleDPIAdjust();
    helpPage();
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

    for (int i = 0; i < fem_system.nodes.size(); ++i)
    {
        ImGui::PushID(i);

        auto &editable_node = fem_system.nodes[i];

        ImGui::Text("Node %d", i + 1);

        // Editable position
        float pos[2] = {editable_node.position[0], editable_node.position[1]};
        if (ImGui::InputFloat2("Position (X, Y) [m]", pos))
        {
            editable_node.position[0] = pos[0];
            editable_node.position[1] = pos[1];
            constraints_changed = true;
        }

        // Constraint type
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
            float angle_deg = editable_node.constraint_angle;
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

            // Remove any beams that reference this node
            for (int s = 0; s < static_cast<int>(fem_system.beams.size()); ++s)
            {
                auto &sp = fem_system.beams[s];
                if (sp.nodes[0] == removed_index || sp.nodes[1] == removed_index)
                {
                    fem_system.beams.erase(fem_system.beams.begin() + s);
                    --s;
                }
            }

            // Decrement node indices in remaining beams that were after the removed node
            for (auto &sp : fem_system.beams)
            {
                if (sp.nodes[0] > removed_index)
                    --sp.nodes[0];
                if (sp.nodes[1] > removed_index)
                    --sp.nodes[1];
            }

            fem_system.nodes.erase(fem_system.nodes.begin() + removed_index);

            ImGui::PopID();
            fem_system.solve_system();
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
            fem_system.nodes.emplace_back(new_x, new_y, Free, 0.0f);
            break;
        case 1:
            fem_system.nodes.emplace_back(new_x, new_y, Fixed, 0.0f);
            break;
        case 2:
            fem_system.nodes.emplace_back(new_x, new_y, FixedPin, 0.0f);
            break;
        case 3:
            fem_system.nodes.emplace_back(new_x, new_y, Slider, new_angle);
            break;
        default:
            fem_system.nodes.emplace_back(new_x, new_y, Free, 0.0f);
            break;
        }

        // Re-solve so system matrices/vectors are rebuilt (solve_system should handle resizing state)
        fem_system.solve_system();
    }

    ImGui::End();

    if (constraints_changed)
    {
        fem_system.solve_system();
    }
}

void GUIHandler::beamEditor()
{
    if (!show_beam_editor)
        return;

    ImGui::Begin("Beam Editor", &show_beam_editor, ImGuiWindowFlags_AlwaysAutoResize);

    // --- Prepare node labels for combos (keep strings alive while UI draws) ---
    std::vector<std::string> node_labels;
    node_labels.reserve(fem_system.nodes.size());
    for (size_t i = 0; i < fem_system.nodes.size(); ++i)
        node_labels.push_back("Node " + std::to_string(i + 1));

    std::vector<const char *> node_items;
    node_items.reserve(node_labels.size());
    for (auto &s : node_labels)
        node_items.push_back(s.c_str());

    // --- Prepare profile & material label arrays ---
    std::vector<const char *> profile_items;
    profile_items.reserve(fem_system.beam_profiles_list.size());
    for (auto &profile : fem_system.beam_profiles_list)
        profile_items.push_back(profile.name.c_str());

    std::vector<const char *> material_items;
    material_items.reserve(fem_system.materials_list.size());
    for (auto &mat : fem_system.materials_list)
        material_items.push_back(mat.name.c_str());

    bool beams_changed = false;

    // --- Existing beams ---
    for (int i = 0; i < static_cast<int>(fem_system.beams.size()); ++i)
    {
        ImGui::PushID(i);
        Beam &beam = fem_system.beams[i];

        ImGui::Text("Beam %d", i + 1);

        // Validate endpoints quickly — if invalid, show warning and allow removal
        const int node_count = static_cast<int>(fem_system.nodes.size());
        bool endpoints_valid = (beam.nodes[0] >= 0 && beam.nodes[0] < node_count &&
                                beam.nodes[1] >= 0 && beam.nodes[1] < node_count);

        if (!endpoints_valid)
        {
            ImGui::TextColored(ImVec4(1.0f, 0.25f, 0.25f, 1.0f), "Invalid endpoints (node index out of range)");
            if (ImGui::Button("Remove Invalid Beam"))
            {
                fem_system.beams.erase(fem_system.beams.begin() + i);
                beams_changed = true;
                ImGui::PopID();
                --i;
                continue;
            }
            ImGui::Separator();
            ImGui::PopID();
            continue;
        }

        // Node A combo
        int node_a = beam.nodes[0];
        if (!node_items.empty())
        {
            if (ImGui::Combo("Node A", &node_a, node_items.data(), static_cast<int>(node_items.size())))
            {
                beam.nodes[0] = node_a;
                beams_changed = true;
            }
        }
        else
        {
            ImGui::TextDisabled("No nodes available");
        }

        // Node B combo
        int node_b = beam.nodes[1];
        if (!node_items.empty())
        {
            if (ImGui::Combo("Node B", &node_b, node_items.data(), static_cast<int>(node_items.size())))
            {
                beam.nodes[1] = node_b;
                beams_changed = true;
            }
        }

        // Prevent Node A == Node B
        if (beam.nodes[0] == beam.nodes[1])
        {
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "Invalid: Node A == Node B");
        }

        // Profile combo (uses index stored in beam.shape_idx)
        int profile_idx = (beam.shape_idx >= 0 && beam.shape_idx < static_cast<int>(fem_system.beam_profiles_list.size()))
                              ? beam.shape_idx
                              : 0;
        if (!profile_items.empty())
        {
            if (ImGui::Combo("Profile", &profile_idx, profile_items.data(), static_cast<int>(profile_items.size())))
            {
                beam.shape_idx = profile_idx;
                beams_changed = true;
            }
        }
        else
        {
            ImGui::TextDisabled("No profiles available");
        }

        // Material combo (uses index stored in beam.material_idx)
        int material_idx = (beam.material_idx >= 0 && beam.material_idx < static_cast<int>(fem_system.materials_list.size()))
                               ? beam.material_idx
                               : 0;
        if (!material_items.empty())
        {
            if (ImGui::Combo("Material", &material_idx, material_items.data(), static_cast<int>(material_items.size())))
            {
                beam.material_idx = material_idx;
                beams_changed = true;
            }
        }
        else
        {
            ImGui::TextDisabled("No materials available");
        }

        // Display stress
        ImGui::SameLine();
        ImGui::Text("Stress: %.2f", beam.stress);

        // Remove button
        if (ImGui::Button("Remove Beam"))
        {
            fem_system.beams.erase(fem_system.beams.begin() + i);
            beams_changed = true;
            ImGui::PopID();
            --i;
            continue;
        }

        ImGui::Separator();
        ImGui::PopID();
    }

    // --- Add new beam controls ---
    static int new_node_a = 0;
    static int new_node_b = 0;
    static int new_profile_idx = 0;
    static int new_material_idx = 0;

    if (!node_items.empty())
    {
        ImGui::Text("Create New Beam:");
        ImGui::Combo("New Node A", &new_node_a, node_items.data(), static_cast<int>(node_items.size()));
        ImGui::Combo("New Node B", &new_node_b, node_items.data(), static_cast<int>(node_items.size()));

        if (!profile_items.empty())
            ImGui::Combo("Profile##NewBeam", &new_profile_idx, profile_items.data(), static_cast<int>(profile_items.size()));
        else
            ImGui::TextDisabled("No profiles available");

        if (!material_items.empty())
            ImGui::Combo("Material##NewBeam", &new_material_idx, material_items.data(), static_cast<int>(material_items.size()));
        else
            ImGui::TextDisabled("No materials available");

        if (new_node_a == new_node_b)
        {
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "Node A and B must be different to create a beam.");
        }

        ImGui::SameLine();
        if (ImGui::Button("Add Beam"))
        {
            // Validate before creating
            if (new_node_a != new_node_b &&
                new_node_a >= 0 && new_node_a < static_cast<int>(fem_system.nodes.size()) &&
                new_node_b >= 0 && new_node_b < static_cast<int>(fem_system.nodes.size()) &&
                !profile_items.empty() && !material_items.empty())
            {
                // Constructor assumed Beam(int n1, int n2, int material_idx, int shape_idx)
                fem_system.beams.emplace_back(new_node_a, new_node_b, new_material_idx, new_profile_idx);
                beams_changed = true;
            }
        }
    }
    else
    {
        ImGui::TextDisabled("No nodes available to create beams.");
    }

    ImGui::End();

    if (beams_changed)
    {
        // Sanity-check beams before solving: remove any that are invalid (defensive)
        for (int i = 0; i < static_cast<int>(fem_system.beams.size()); ++i)
        {
            const Beam &s = fem_system.beams[i];
            if (s.nodes[0] < 0 || s.nodes[0] >= static_cast<int>(fem_system.nodes.size()) ||
                s.nodes[1] < 0 || s.nodes[1] >= static_cast<int>(fem_system.nodes.size()) ||
                s.nodes[0] == s.nodes[1] ||
                s.material_idx < 0 || s.material_idx >= static_cast<int>(fem_system.materials_list.size()) ||
                s.shape_idx < 0 || s.shape_idx >= static_cast<int>(fem_system.beam_profiles_list.size()))
            {
                std::cerr << "Removing invalid beam at index " << i << " during validation.\n";
                fem_system.beams.erase(fem_system.beams.begin() + i);
                --i;
            }
        }

        fem_system.solve_system();
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
    for (int i = 0; i < fem_system.nodes.size(); ++i)
    {
        if (fem_system.nodes[i].constraint_type == Free || fem_system.nodes[i].constraint_type == Slider)
        {
            float fx = fem_system.forces(i * 3);
            float fy = fem_system.forces(i * 3 + 1);
            ImGui::PushID(i);

            ImGui::Text("Node %d Forces:", i + 1);

            // Fx controls
            ImGui::Text("Fx:");
            ImGui::SameLine();
            if (ImGui::SliderFloat("##FxSlider", &fx, -10000.0f, 10000.0f))
            {
                fem_system.forces(i * 3) = fx;
                forces_changed = true;
            }
            ImGui::SameLine();
            if (ImGui::InputFloat("##FxInput", &fx))
            {
                fem_system.forces(i * 3) = fx;
                forces_changed = true;
            }

            // Fy controls
            ImGui::Text("Fy:");
            ImGui::SameLine();
            if (ImGui::SliderFloat("##FySlider", &fy, -10000.0f, 10000.0f))
            {
                fem_system.forces(i * 3 + 1) = fy;
                forces_changed = true;
            }
            ImGui::SameLine();
            if (ImGui::InputFloat("##FyInput", &fy))
            {
                fem_system.forces(i * 3 + 1) = fy;
                forces_changed = true;
            }

            ImGui::Separator();
            ImGui::PopID();
        }
    }

    if (forces_changed)
        fem_system.solve_system();

    // Solution display
    ImGui::Text("Solution:");
    for (int i = 0; i < fem_system.nodes.size(); ++i)
    {
        // Convert rotation from radians to degrees for display
        double theta_rad = fem_system.displacement(i * 3 + 2);
        double theta_deg = theta_rad * 180.0 / M_PI; // Assuming M_PI is defined

        ImGui::Text("Node %d: u=%.6f m, v=%.6f m, theta=%.6f deg", i + 1,
                    fem_system.displacement(i * 3),     // NEW INDEX
                    fem_system.displacement(i * 3 + 1), // NEW INDEX
                    theta_deg);                         // NEW: Rotation
    }

    // Beam stresses
    ImGui::Text("Beam Stresses (MPa):");
    int index = 1;
    for (const auto &beam : fem_system.beams)
    {
        sf::Color color = renderer.getStressColor(beam.stress, fem_system.min_stress, fem_system.max_stress);
        ImGui::TextColored(ImVec4(color.r / 255.f, color.g / 255.f, color.b / 255.f, 1.f),
                           "Beam %d: %.2f MPa", index++, beam.stress);
    }

    ImGui::End();
}

void GUIHandler::handleLoadPopup()
{
    if (request_load_popup)
    {
        ImGui::OpenPopup("Load From");
        request_load_popup = false;
    }

    if (ImGui::BeginPopupModal("Load From", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
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
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    // Display error message if load failed
    if (load_error)
    {
        ImGui::OpenPopup("Load Error");
    }
    if (ImGui::BeginPopupModal("Load Error", &load_error, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Failed to load file: %s", error_msg.c_str());
        if (ImGui::Button("OK"))
        {
            load_error = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    if (!trigger_load_read)
        return;

    trigger_load_read = false;
    error_msg.clear();

    // check if user put .ffem extension, add if missing
    std::string filename_str(filename_buf);
    if (filename_str.size() < 5 || filename_str.substr(filename_str.size() - 5) != ".ffem")
    {
        filename_str += ".ffem";
        std::strncpy(filename_buf, filename_str.c_str(), sizeof(filename_buf));
        filename_buf[sizeof(filename_buf) - 1] = '\0'; // ensure null termination
    }

    std::ifstream ifs(filename_buf, std::ios::binary);
    if (!ifs)
    {
        error_msg = "Could not open file for reading.";
        load_error = true;
        return;
    }

    // 1. Header/Version
    uint32_t version = 0;
    ifs.read(reinterpret_cast<char *>(&version), sizeof(version));
    if (!ifs || version != FILE_VERSION)
    {
        error_msg = "File version mismatch or read error. Expected " + std::to_string(FILE_VERSION) + ", got " + std::to_string(version);
        load_error = true;
        return;
    }

    // Clear current system
    fem_system.materials_list.clear();
    fem_system.beam_profiles_list.clear();
    fem_system.nodes.clear();
    fem_system.beams.clear();
    fem_system.forces = Eigen::VectorXd();

    // 2. Material Profiles
    uint32_t material_count = 0;
    ifs.read(reinterpret_cast<char *>(&material_count), sizeof(material_count));
    if (!ifs)
    {
        error_msg = "Failed reading material count.";
        load_error = true;
        return;
    }
    fem_system.materials_list.resize(material_count);
    for (uint32_t i = 0; i < material_count; ++i)
    {
        MaterialProfile &m = fem_system.materials_list[i];
        m.name = readString(ifs);
        ifs.read(reinterpret_cast<char *>(&m.youngs_modulus), sizeof(m.youngs_modulus));
        if (!ifs)
        {
            error_msg = "Failed reading material profile " + std::to_string(i);
            load_error = true;
            return;
        }
    }

    // 3. Beam Profiles
    uint32_t profile_count = 0;
    ifs.read(reinterpret_cast<char *>(&profile_count), sizeof(profile_count));
    if (!ifs)
    {
        error_msg = "Failed reading beam profile count.";
        load_error = true;
        return;
    }
    fem_system.beam_profiles_list.resize(profile_count);
    for (uint32_t i = 0; i < profile_count; ++i)
    {
        BeamProfile &p = fem_system.beam_profiles_list[i];
        p.name = readString(ifs);
        ifs.read(reinterpret_cast<char *>(&p.area), sizeof(p.area));
        ifs.read(reinterpret_cast<char *>(&p.moment_of_inertia), sizeof(p.moment_of_inertia));
        ifs.read(reinterpret_cast<char *>(&p.section_modulus), sizeof(p.section_modulus));
        if (!ifs)
        {
            error_msg = "Failed reading beam profile " + std::to_string(i);
            load_error = true;
            return;
        }
    }

    // 4. Nodes
    uint32_t node_count = 0;
    ifs.read(reinterpret_cast<char *>(&node_count), sizeof(node_count));
    if (!ifs)
    {
        error_msg = "Failed reading node count.";
        load_error = true;
        return;
    }
    fem_system.nodes.resize(node_count);
    for (uint32_t i = 0; i < node_count; ++i)
    {
        Node &n = fem_system.nodes[i];
        ifs.read(reinterpret_cast<char *>(&n.position[0]), sizeof(float) * 2);
        int32_t type = 0;
        ifs.read(reinterpret_cast<char *>(&type), sizeof(type));
        ifs.read(reinterpret_cast<char *>(&n.constraint_angle), sizeof(n.constraint_angle));
        if (!ifs)
        {
            error_msg = "Failed reading node " + std::to_string(i);
            load_error = true;
            return;
        }
        // validate enum range (defensive)
        if (type < static_cast<int32_t>(Free) || type > static_cast<int32_t>(Slider))
            n.constraint_type = Free;
        else
            n.constraint_type = static_cast<ConstraintType>(type);
    }

    // 5. Beams (now using indices)
    uint32_t beam_count = 0;
    ifs.read(reinterpret_cast<char *>(&beam_count), sizeof(beam_count));
    if (!ifs)
    {
        error_msg = "Failed reading beam count.";
        load_error = true;
        return;
    }
    fem_system.beams.resize(beam_count);
    for (uint32_t i = 0; i < beam_count; ++i)
    {
        Beam &s = fem_system.beams[i];

        // Nodes indices (int32_t[2])
        int32_t n0 = -1, n1 = -1;
        ifs.read(reinterpret_cast<char *>(&n0), sizeof(n0));
        ifs.read(reinterpret_cast<char *>(&n1), sizeof(n1));
        if (!ifs)
        {
            error_msg = "Failed reading beam node indices for beam " + std::to_string(i);
            load_error = true;
            return;
        }
        s.nodes[0] = n0;
        s.nodes[1] = n1;

        // Stress (float)
        ifs.read(reinterpret_cast<char *>(&s.stress), sizeof(s.stress));
        if (!ifs)
        {
            error_msg = "Failed reading beam stress for beam " + std::to_string(i);
            load_error = true;
            return;
        }

        // Material and Profile Indices (int32_t[2])
        int32_t mat_idx = -1, shape_idx = -1;
        ifs.read(reinterpret_cast<char *>(&mat_idx), sizeof(mat_idx));
        ifs.read(reinterpret_cast<char *>(&shape_idx), sizeof(shape_idx));
        if (!ifs)
        {
            error_msg = "Failed reading beam material/shape indices for beam " + std::to_string(i);
            load_error = true;
            return;
        }

        // Validate indices
        if (mat_idx < 0 || mat_idx >= static_cast<int32_t>(fem_system.materials_list.size()) ||
            shape_idx < 0 || shape_idx >= static_cast<int32_t>(fem_system.beam_profiles_list.size()))
        {
            error_msg = "Invalid material/shape index in beam " + std::to_string(i);
            load_error = true;
            return;
        }
        s.material_idx = mat_idx;
        s.shape_idx = shape_idx;

        // validate node indices
        if (s.nodes[0] < 0 || s.nodes[0] >= static_cast<int32_t>(fem_system.nodes.size()) ||
            s.nodes[1] < 0 || s.nodes[1] >= static_cast<int32_t>(fem_system.nodes.size()))
        {
            error_msg = "Invalid node index in beam " + std::to_string(i);
            load_error = true;
            return;
        }
    }

    // 6. Forces (Eigen::VectorXd)
    uint32_t fcount = 0;
    ifs.read(reinterpret_cast<char *>(&fcount), sizeof(fcount));
    if (!ifs)
    {
        error_msg = "Failed reading forces count.";
        load_error = true;
        return;
    }

    const size_t expected_forces = fem_system.nodes.size() * 3;
    if (fcount == expected_forces)
    {
        fem_system.forces.resize(fcount);
        if (fcount > 0)
        {
            ifs.read(reinterpret_cast<char *>(fem_system.forces.data()), sizeof(double) * fcount);
            if (!ifs)
            {
                error_msg = "Failed reading full forces data.";
                load_error = true;
                return;
            }
        }
    }
    else
    {
        // mismatch: safer to zero forces rather than trust a mismatched file
        std::cerr << "Warning: forces count in file (" << fcount << ") does not match expected (" << expected_forces << "). Zeroing forces.\n";
        fem_system.forces = Eigen::VectorXd::Zero(expected_forces);

        // Attempt to skip over the saved forces blob if file contains it
        if (fcount > 0)
        {
            // seek forward to skip the raw doubles (best-effort)
            ifs.seekg(static_cast<std::streamoff>(sizeof(double) * fcount), std::ios::cur);
            // ignore seek failure — we already zeroed forces
        }
    }

    // final file-read sanity check
    if (ifs.fail() && !ifs.eof())
    {
        error_msg = "Error occurred during file reading or file truncated.";
        load_error = true;
        return;
    }

    fem_system.total_dof = static_cast<int>(fem_system.nodes.size()) * 3;
    fem_system.displacement = Eigen::VectorXd::Zero(fem_system.total_dof);
    // forces already set (either loaded or zeroed)
    fem_system.solve_system();
    renderer.autoZoomToFit();

    // all done
}

void GUIHandler::handleSavePopup()
{
    if (request_save_popup)
    {
        ImGui::OpenPopup("Save As");
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
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    }

    // Display error message if save failed
    if (save_error)
    {
        ImGui::OpenPopup("Save Error");
    }
    if (ImGui::BeginPopupModal("Save Error", &save_error, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Failed to save file: %s", error_msg.c_str());
        if (ImGui::Button("OK"))
        {
            save_error = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    if (!trigger_save_write)
        return;

    trigger_save_write = false;
    error_msg.clear();

    // check if user put .ffem extension, add if missing
    std::string filename_str(filename_buf);
    if (filename_str.size() < 5 || filename_str.substr(filename_str.size() - 5) != ".ffem")
    {
        filename_str += ".ffem";
        std::strncpy(filename_buf, filename_str.c_str(), sizeof(filename_buf));
        filename_buf[sizeof(filename_buf) - 1] = '\0'; // ensure null termination
    }

    std::ofstream ofs(filename_buf, std::ios::binary);
    if (!ofs)
    {
        error_msg = "Could not open file for writing.";
        save_error = true;
        return;
    }

    // 1. Header/Version
    ofs.write(reinterpret_cast<const char *>(&FILE_VERSION), sizeof(FILE_VERSION));

    // 2. Material Profiles
    uint32_t material_count = static_cast<uint32_t>(fem_system.materials_list.size());
    ofs.write(reinterpret_cast<const char *>(&material_count), sizeof(material_count));
    for (const auto &m : fem_system.materials_list)
    {
        writeString(ofs, m.name);
        ofs.write(reinterpret_cast<const char *>(&m.youngs_modulus), sizeof(m.youngs_modulus));
    }

    // 3. Beam Profiles
    uint32_t profile_count = static_cast<uint32_t>(fem_system.beam_profiles_list.size());
    ofs.write(reinterpret_cast<const char *>(&profile_count), sizeof(profile_count));
    for (const auto &p : fem_system.beam_profiles_list)
    {
        writeString(ofs, p.name);
        ofs.write(reinterpret_cast<const char *>(&p.area), sizeof(p.area));
        ofs.write(reinterpret_cast<const char *>(&p.moment_of_inertia), sizeof(p.moment_of_inertia));
        ofs.write(reinterpret_cast<const char *>(&p.section_modulus), sizeof(p.section_modulus));
    }

    // 4. Nodes
    uint32_t node_count = static_cast<uint32_t>(fem_system.nodes.size());
    ofs.write(reinterpret_cast<const char *>(&node_count), sizeof(node_count));
    for (const auto &n : fem_system.nodes)
    {
        ofs.write(reinterpret_cast<const char *>(&n.position[0]), sizeof(float) * 2);
        int32_t type = static_cast<int32_t>(n.constraint_type);
        ofs.write(reinterpret_cast<const char *>(&type), sizeof(type));
        ofs.write(reinterpret_cast<const char *>(&n.constraint_angle), sizeof(n.constraint_angle));
    }

    // 5. Beams (write node indices, stress, material_idx, shape_idx)
    uint32_t beam_count = static_cast<uint32_t>(fem_system.beams.size());
    ofs.write(reinterpret_cast<const char *>(&beam_count), sizeof(beam_count));
    for (const auto &s : fem_system.beams)
    {
        // Nodes indices (int32_t[2])
        int32_t n0 = s.nodes[0];
        int32_t n1 = s.nodes[1];
        ofs.write(reinterpret_cast<const char *>(&n0), sizeof(n0));
        ofs.write(reinterpret_cast<const char *>(&n1), sizeof(n1));

        // Stress (float)
        ofs.write(reinterpret_cast<const char *>(&s.stress), sizeof(s.stress));

        // Material and Profile Indices (int32_t[2])
        int32_t mat_idx = s.material_idx;
        int32_t shape_idx = s.shape_idx;
        ofs.write(reinterpret_cast<const char *>(&mat_idx), sizeof(mat_idx));
        ofs.write(reinterpret_cast<const char *>(&shape_idx), sizeof(shape_idx));
    }

    // 6. Forces (Eigen::VectorXd)
    uint32_t fcount = static_cast<uint32_t>(fem_system.forces.size());
    ofs.write(reinterpret_cast<const char *>(&fcount), sizeof(fcount));
    if (fcount > 0)
    {
        ofs.write(reinterpret_cast<const char *>(fem_system.forces.data()), sizeof(double) * fcount);
    }

    if (ofs.fail())
    {
        error_msg = "Error occurred during file writing.";
        save_error = true;
        return;
    }

    // success - file closed on destruction of ofs
}

void GUIHandler::profileEditor()
{
    if (!show_profile_editor)
        return;

    ImGui::Begin("Profile Editor", &show_profile_editor, ImGuiWindowFlags_AlwaysAutoResize);

    bool profiles_changed = false;

    // --- Existing profiles ---
    if (fem_system.beam_profiles_list.empty())
    {
        ImGui::TextDisabled("No profiles created.");
    }
    else
    {
        for (int i = 0; i < static_cast<int>(fem_system.beam_profiles_list.size()); ++i)
        {
            ImGui::PushID(i);
            BeamProfile &profile = fem_system.beam_profiles_list[i];

            ImGui::Text("Profile %d: %s", i + 1, profile.name.c_str());

            // --- Editable Area ---
            float area = profile.area;
            if (ImGui::InputFloat("Area (m^2)", &area))
            {
                profile.area = area;

                for (auto &b : fem_system.beams)
                {
                    if (b.shape_idx == i)
                        fem_system.beam_profiles_list[b.shape_idx].area = area;
                }
                profiles_changed = true;
            }

            // --- Editable Moment of Inertia ---
            float I = profile.moment_of_inertia;
            if (ImGui::InputFloat("Moment of Inertia I (m^4)", &I))
            {
                profile.moment_of_inertia = I;

                for (auto &b : fem_system.beams)
                {
                    if (b.shape_idx == i)
                        fem_system.beam_profiles_list[b.shape_idx].moment_of_inertia = I;
                }
                profiles_changed = true;
            }

            // --- Editable Section Modulus ---
            float S = profile.section_modulus;
            if (ImGui::InputFloat("Section Modulus S (m^3)", &S))
            {
                profile.section_modulus = S;

                for (auto &b : fem_system.beams)
                {
                    if (b.shape_idx == i)
                        fem_system.beam_profiles_list[b.shape_idx].section_modulus = S;
                }
                profiles_changed = true;
            }

            // Remove profile button
            ImGui::SameLine();
            if (ImGui::Button("Remove Profile"))
            {
                // Remove beams using this profile
                for (int s = 0; s < static_cast<int>(fem_system.beams.size()); ++s)
                {
                    if (fem_system.beams[s].shape_idx == i)
                    {
                        fem_system.beams.erase(fem_system.beams.begin() + s);
                        --s;
                    }
                }

                fem_system.beam_profiles_list.erase(fem_system.beam_profiles_list.begin() + i);

                ImGui::PopID();
                profiles_changed = true;
                --i;
                continue;
            }

            ImGui::Separator();
            ImGui::PopID();
        }
    }

    // --- New Profile Creation ---
    ImGui::Separator();
    ImGui::Text("Create New Profile:");

    static char new_profile_name[128] = "";
    static float new_profile_area = 0.1963f;
    static float new_profile_I = 0.005f;
    static float new_profile_S = 0.01f;

    ImGui::InputText("Name", new_profile_name, sizeof(new_profile_name));
    ImGui::InputFloat("Area (m^2)", &new_profile_area);
    ImGui::InputFloat("Moment of Inertia I (m^4)", &new_profile_I);
    ImGui::InputFloat("Section Modulus S (m^3)", &new_profile_S);

    ImGui::SameLine();
    if (ImGui::Button("Add Profile"))
    {
        std::string name_str(new_profile_name);
        bool duplicate = false;

        for (const auto &p : fem_system.beam_profiles_list)
        {
            if (p.name == name_str)
            {
                duplicate = true;
                break;
            }
        }

        if (!name_str.empty() && !duplicate)
        {
            BeamProfile new_p;
            new_p.name = name_str;
            new_p.area = new_profile_area;
            new_p.moment_of_inertia = new_profile_I;
            new_p.section_modulus = new_profile_S;

            fem_system.beam_profiles_list.push_back(new_p);

            // Reset fields
            new_profile_name[0] = '\0';
            new_profile_area = 0.1963f;
            new_profile_I = 0.005f;
            new_profile_S = 0.01f;

            profiles_changed = true;
        }
    }

    ImGui::End();

    if (profiles_changed)
    {
        fem_system.solve_system();
    }
}

void GUIHandler::materialEditor()
{
    if (!show_material_editor)
        return;

    ImGui::Begin("Material Editor", &show_material_editor, ImGuiWindowFlags_AlwaysAutoResize);

    bool materials_changed = false;

    // --- Existing materials ---
    if (fem_system.materials_list.empty())
    {
        ImGui::TextDisabled("No materials created.");
    }
    else
    {
        for (int i = 0; i < static_cast<int>(fem_system.materials_list.size()); ++i)
        {
            ImGui::PushID(i);
            MaterialProfile &mat = fem_system.materials_list[i];

            ImGui::Text("Material %d: %s", i + 1, mat.name.c_str());

            // Editable Young's modulus
            float youngs = mat.youngs_modulus;
            if (ImGui::InputFloat("Young's Modulus (Pa)", &youngs))
            {
                mat.youngs_modulus = youngs;
                // Update all beams using this material index
                for (auto &sp : fem_system.beams)
                {
                    if (sp.material_idx == i)
                        sp.material_idx = i; // Index remains the same, modulus updated via material reference
                }
                materials_changed = true;
            }

            // Remove material button
            ImGui::SameLine();
            if (ImGui::Button("Remove Material"))
            {
                // Remove all beams using this material
                for (int s = 0; s < static_cast<int>(fem_system.beams.size()); ++s)
                {
                    if (fem_system.beams[s].material_idx == i)
                    {
                        fem_system.beams.erase(fem_system.beams.begin() + s);
                        --s;
                    }
                }

                // Erase the material
                fem_system.materials_list.erase(fem_system.materials_list.begin() + i);

                ImGui::PopID();
                materials_changed = true;

                --i; // Adjust loop index after erase
                continue;
            }

            ImGui::Separator();
            ImGui::PopID();
        }
    }

    // --- New material creation ---
    ImGui::Separator();
    ImGui::Text("Create New Material:");
    static char new_mat_name[128] = "";
    static float new_mat_youngs = 30e6f;

    ImGui::InputText("Name", new_mat_name, sizeof(new_mat_name));
    ImGui::InputFloat("Young's Modulus (Pa)", &new_mat_youngs);

    ImGui::SameLine();
    if (ImGui::Button("Add Material"))
    {
        std::string name_str(new_mat_name);
        bool duplicate = false;
        for (const auto &mat : fem_system.materials_list)
        {
            if (mat.name == name_str)
            {
                duplicate = true;
                break;
            }
        }

        if (!name_str.empty() && !duplicate)
        {
            MaterialProfile new_mat;
            new_mat.name = name_str;
            new_mat.youngs_modulus = new_mat_youngs;
            fem_system.materials_list.push_back(new_mat);

            // Reset fields
            new_mat_name[0] = '\0';
            new_mat_youngs = 30e6f;

            materials_changed = true;
        }
    }

    ImGui::End();

    if (materials_changed)
    {
        fem_system.solve_system();
    }
}

void GUIHandler::handleDPIAdjust()
{
    static bool popup_just_opened = false;

    if (request_dpi_adjust)
    {
        request_dpi_adjust = false;
        popup_just_opened = true; // mark that we need to init on next BeginPopupModal
        ImGui::OpenPopup("Adjust DPI Scaling");
    }

    static float dpi_percent = current_dpi_scale * 100.0f;

    if (ImGui::BeginPopupModal("Adjust DPI Scaling", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (popup_just_opened)
        {
            // Initialize only ONCE when popup opens
            float detected = current_dpi_scale;

            try
            {
                detected = renderer.GetDPIScale(window.getNativeHandle());
            }
            catch (...)
            {
            }

            dpi_percent = detected * 100.0f;
            popup_just_opened = false;
        }

        // Close button
        float close_btn_sz = ImGui::GetFrameHeight();
        ImGui::SameLine(ImGui::GetWindowWidth() - close_btn_sz - ImGui::GetStyle().FramePadding.x);
        if (ImGui::Button("X"))
            ImGui::CloseCurrentPopup();

        float detected = current_dpi_scale;
        try
        {
            detected = renderer.GetDPIScale(window.getNativeHandle());
        }
        catch (...)
        {
        }

        ImGui::Text("Detected DPI Scale: %.2f (%.0f%%)", detected, detected * 100.0f);
        ImGui::SliderFloat("Scale (%)", &dpi_percent, 10.0f, 400.0f, "%.0f%%");
        ImGui::TextWrapped("Apply changes to scale ImGui fonts and widgets. 100% = no scaling.");

        if (ImGui::Button("Apply"))
        {
            // DON'T call applyDPIScale here.
            // Just store the request for the start of the next frame.
            pending_dpi_scale = dpi_percent / 100.0f;

            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel"))
            ImGui::CloseCurrentPopup();

        ImGui::EndPopup();
    }
}

void GUIHandler::headerBar()
{
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("New System", "Ctrl+N"))
            {
                fem_system.nodes.clear();
                fem_system.beams.clear();
                fem_system.solve_system();
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
            if (ImGui::MenuItem("Adjust DPI Scaling"))
            {
                request_dpi_adjust = true;
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
            if (ImGui::MenuItem("System Beam Editor"))
            {
                show_beam_editor = !show_beam_editor;
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

            ImGui::EndMenu();
        }
        if (ImGui::MenuItem("Help"))
        {
            show_help_page = !show_help_page;
        }

        {
            float frame_rate = ImGui::GetIO().Framerate;
            char fps_buf[32];
            std::snprintf(fps_buf, sizeof(fps_buf), "FPS: %.1f", frame_rate);

            // Compute text width and place it so it never goes past the right window padding
            float text_w = ImGui::CalcTextSize(fps_buf).x;
            ImGuiStyle &style = ImGui::GetStyle();
            float window_w = ImGui::GetWindowWidth();
            float target_x = window_w - text_w - style.WindowPadding.x - style.ItemSpacing.x;

            // Ensure we don't move left of the current cursor position
            float cur_x = ImGui::GetCursorPosX();
            if (target_x < cur_x)
                target_x = cur_x;

            ImGui::SetCursorPosX(target_x);
            ImGui::Text("%s", fps_buf);
        }

        ImGui::EndMainMenuBar();
    }
}

void GUIHandler::helpPage()
{
    if (!show_help_page)
        return;
    ImGui::Begin("Help", &show_help_page, ImGuiWindowFlags_AlwaysAutoResize);

    ImGui::Text("2D Truss System Solver - Help");
    ImGui::Separator();

    if (ImGui::CollapsingHeader("About", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::TextWrapped("This application allows you to create and analyze 2D truss systems using the finite element method. You can add nodes and beams, define material properties and beam profiles, apply forces, and view the resulting displacements and stresses.");
    }

    if (ImGui::CollapsingHeader("How to Use", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::TextWrapped("1. Use the 'System Node Editor' to add and manage nodes in your truss system. You can set constraints for each node (fixed, free, or slider).");
        ImGui::TextWrapped("2. Use the 'System Beam Editor' to add beams between nodes. Select material profiles and beam profiles for each beam.");
        ImGui::TextWrapped("3. Use the 'Material Editor' to create and manage material profiles with specific Young's modulus values.");
        ImGui::TextWrapped("4. Use the 'Profile Editor' to create and manage beam profiles with specific cross-sectional areas.");
        ImGui::TextWrapped("5. Use the 'System Controls' to apply forces to nodes and view the computed displacements and stresses.");
    }

    if (ImGui::CollapsingHeader("Extra Features", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // Extra features as bullet list
        ImGui::BulletText("Save and load truss systems using the 'File' menu. These use the custom .ffem binary format.");
        ImGui::BulletText("Adjust DPI scaling for better visibility on high-resolution displays.");
        ImGui::BulletText("Auto view and zoom features to fit the truss system within the viewport.");
    }

    ImGui::End();
}