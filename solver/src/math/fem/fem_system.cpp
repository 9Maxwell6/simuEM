#include "math/fem/fem_system.h"

using namespace simu;



FEM_System::FEM_System(Mesh& mesh):mesh_(mesh){
    dim_ = mesh_.get_mesh_dimension();

};



std::vector<FEM_Space *> FEM_System::create_FE_space(Space fs, int p_order){

    std::vector<FEM_Space *> new_fem_space;

    switch (dim_)
    {
    case 3:
        if(mesh_.)
        break;
    case 2:
        break;
    case 1:
        break;
    default:
        Logger::error("FEM_System::create_FE_space - impossible dimension: "+std::to_string(dim_));
        
    }
}


/**
 * @brief Assign space to all elements in mesh
 *
 * @param fs finite element space.
 * @param p_order order of basis function.
 */
void FEM_System::assign_FE_space(Space fs, int p_order)
{
    // TODO: should use vector<Space>, what if we assign duplicate Space?
    // use vector<FEM_Space *>,
    //global_space.push_back(fs);
    
}
