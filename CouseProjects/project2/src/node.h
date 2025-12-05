#pragma once

enum ConstraintType
{
    Free,
    Fixed,
    FixedPin,
    Slider
};

class Node
{
public:
    float position[2];
    ConstraintType constraint_type;
    float constraint_angle;

    // Default constructor required for std::vector::resize
    Node()
        : position{0.0f, 0.0f},
          constraint_type(Free),
          constraint_angle(0.0f)
    {
    }

    Node(float x, float y, ConstraintType ct = Free, float angle = 0.0f)
        : position{x, y},
          constraint_type(ct),
          constraint_angle(angle)
    {
    }
};