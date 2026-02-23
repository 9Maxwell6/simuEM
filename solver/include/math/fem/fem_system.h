#pragma once

#include "entity/mesh/mesh.h"
#include "math/fem/space_collection.h"
#include "entity/mesh/e_collection.h"
#include "math/fem/dof_handler.h"



namespace simu {

class FEM_System
{

private:
    int dim_;

    Mesh& mesh_;
    DoF_Handler dof_handler_;

    std::vector<size_t> dof_list_;     // entry -> index in global dof

    std::vector<size_t> dof_space_offset_;   // used for block matrix assemble.

    std::vector<size_t> elements_dof_lookup_list;  // [dofs of element 1, dofs of element 2, ...]

    std::vector<FEM_Space * > global_space;
    std::unordered_map<Key, std::vector<FEM_Space * >, Key_Hash> group_space;
    



    //std::vector<size_t> elements_dof_lookup_list;

    void assign_dof(Element * e);

    bool initialize_FE_space(FEM_Space& fe_space);

    static Basis_Shape to_basis_shape(Geometry t);

    static Geometry to_element_geometry(Basis_Shape t);


public:
    FEM_System(Mesh& mesh);

    // assign functional space to all elements in mesh.
    bool assign_FE_space(FEM_Space& fe_space);

    // assign functional space to specific group of elements.
    bool assign_group_FE_space(FEM_Space& fe_space, const Key group_key={0,0});




};

}