#pragma once

#include "world/mesh/group_key.h"

#include <vector>

namespace simu {


struct Boundary 
{
    int id;
    Key* b_group;

    Boundary* base_boundary; 

    std::vector<Region*>   region_list;
    std::vector<Boundary*> sub_boundary_list;
};


struct Region 
{
    int id;
    Key* r_group;

    Region* base_region;  

    std::vector<Boundary*> boundary_list;
    std::vector<region*>   sub_region_list;
    
};

}