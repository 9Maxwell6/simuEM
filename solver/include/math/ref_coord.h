#pragma once

struct Ref_Coord
{
    double x, y, z;
};


// integration point on reference element
struct Integration_Point 
{
    Ref_Coord coord;
    double weight;
};