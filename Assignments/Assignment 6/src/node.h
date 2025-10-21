#pragma once

enum ConstraintType
{
    Free,
    Fixed,
    Slider
};

class Node
{
public:
    int id;
    float position[2]; // position in 2D space (x, y)
    ConstraintType constraint_type;
    float constraint_angle; // Angle in degrees (0 for horizontal, 90 for vertical)

    Node(int i, float x, float y, ConstraintType ct = Free, float angle = 0.0f) : id(i), constraint_type(ct), constraint_angle(angle)
    {
        position[0] = x;
        position[1] = y;
    }
};