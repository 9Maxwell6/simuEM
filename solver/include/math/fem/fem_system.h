#pragma once

#include "entity/mesh/mesh.h"
#include "math/fem/space_collection.h"
#include "entity/mesh/e_collection.h"
#include "math/fem/dof_handler.h"



namespace simu {

class FEM_System
{

private:
    Mesh& mesh_;
    DoF_Handler dof_handler_;

    std::vector<size_t> dof_list_;     // entry -> index in global dof

    std::vector<size_t> dof_space_offset_;

    std::vector<size_t> elements_dof_lookup_list;  // [dofs of element 1, dofs of element 2, ...]

    std::vector<FEM_Space *> global_space;
    std::unordered_map<Key, FEM_Space *, Key_Hash> group_space;


    //std::vector<size_t> elements_dof_lookup_list;

    void assign_dof(Element * e);

public:
    FEM_System(Mesh& mesh);

    // assign functional space to all elements in mesh.
    void assign_FE_space(Space fs, int p_order=1);

    // assign functional space to specific group of elements.
    void assign_group_FE_space(Space fs, const Key group_key={0,0}, int p_order=1);




};

}