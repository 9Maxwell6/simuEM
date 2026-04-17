#pragma once

#include "world/mesh/group_key.h"

#include <vector>

namespace simu {

struct Region;

struct Boundary 
{
    int id;
    std::string description;

    Key b_group;

    const Boundary* base_boundary; 

    std::vector<const Region*>   region_list;
    std::vector<const Boundary*> sub_boundary_list;
};


struct Region 
{
    int id;
    std::string description;

    Key r_group;

    const Region* base_region;  

    std::vector<const Boundary*> boundary_list;
    std::vector<const Region*>   sub_region_list;
    
};

}