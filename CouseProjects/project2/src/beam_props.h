#pragma once

struct MaterialProfile
{
    std::string name;
    double youngs_modulus;
};

struct BeamProfile
{
    std::string name;
    double area;
    double moment_of_inertia;
    double section_modulus;
};
