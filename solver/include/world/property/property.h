#pragma once

#include "coefficient.h"

#include <bitset>
#include <stddef.h>
#include <vector>
#include <map>

namespace simu {

struct Property 
{
    std::string name;
    std::vector<double> scalar_data;

    std::vector<double> vector_data_x; 
    std::vector<double> vector_data_y; 
    std::vector<double> vector_data_z; 

    std::vector<Coefficient *> coeff;  
};





class Property_handler
{
private:
    


    //
};

}