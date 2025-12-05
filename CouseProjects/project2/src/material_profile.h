#pragma once

struct MaterialProfile
{
    std::string name;
    double youngs_modulus;
};

class BeamProperties
{
public:
    std::string name;
    double area;
    MaterialProfile material;
};
