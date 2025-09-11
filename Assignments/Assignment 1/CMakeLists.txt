cmake_minimum_required(VERSION 3.16)
project(ImGuiSFMLExample)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# --- Fetch SFML 3.0.0 ---
FetchContent_Declare(
    SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.1
)
FetchContent_MakeAvailable(SFML)

# --- Fetch ImGui (use version required by ImGui-SFML) ---
FetchContent_Declare(
    ImGui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.9  # Use the minimum version required by ImGui-SFML
)
FetchContent_MakeAvailable(ImGui)

# --- Configure ImGui-SFML ---
set(IMGUI_DIR ${imgui_SOURCE_DIR})
set(IMGUI_SFML_FIND_SFML OFF)

# --- Fetch ImGui-SFML (stable version compatible with SFML 2.6.x) ---
FetchContent_Declare(
    ImGui-SFML
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml.git
    GIT_TAG master
)
FetchContent_MakeAvailable(ImGui-SFML)

# --- Your project executable ---
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE 
    ImGui-SFML::ImGui-SFML
    sfml-graphics
    sfml-window
    sfml-system
)