#include "math/fem/fem_system.h"

using namespace simu;



FEM_System::FEM_System(Mesh& mesh):mesh_(mesh)
{
    dim_ = mesh_.get_mesh_dimension();
    dof_offset_ = 0;
};




bool FEM_System::initialize_space_dof()
{

}



size_t FEM_System::assign_element_dof(FEM_Space& fe_space, Element& e)
{
    size_t current_offset = dof_offset_;
    FEM_Space * basis_space = fe_space.get_basis_space(to_basis_shape(e.get_geometry()));
    

    return current_offset;
}




Basis_Shape FEM_System::to_basis_shape(Geometry t)
{
    switch (t) {
        case Geometry::TETRAHEDRON: return Basis_Shape::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - type not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

Geometry FEM_System::to_element_geometry(Basis_Shape g)
{
    switch (g) {
        case Basis_Shape::TETRAHEDRON: return Geometry::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - geometry not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}


/*
bool FEM_System::create_FE_space(FEM_Space * fe_space){

    std::vector<FEM_Space *> new_fem_space;

    switch (dim_)
    {
    case 3:
        //if(mesh_.)
        break;
    case 2:
        break;
    case 1:
        break;
    default:
        Logger::error("FEM_System::create_FE_space - impossible dimension: "+std::to_string(dim_));
        
    }
}
    */


/**
 * @brief Assign space to all elements in mesh
 *
 * @param fs finite element space.
 * @param p_order order of basis function.
 */
//bool FEM_System::initialize_FE_space(FEM_Space& fe_space)
//{
    // TODO: should use vector<Space>, what if we assign duplicate Space?
    // use vector<FEM_Space *>,
    //global_space.push_back(fs);
    
//}


/**
 * @brief Assign space to all elements in mesh
 *
 * @param fs finite element space.
 */
bool FEM_System::assign_FE_space(FEM_Space& fe_space)
{
    for(const auto& [type, size] : mesh_.get_mesh_element_geometry_size())
    {
        fe_space.add_basis_shape(to_basis_shape(type));
    }
    
    // initialize dof
    for(Element* e : mesh_.get_mesh_elements())
    {
        
    }
}