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


    void assign_dof(Element * e);

public:
    FEM_System(Mesh& mesh);

    // assign functional space to all elements in mesh.
    void construct_dof(Space fs, int p_order=1);

    // assign functional space to specific group of elements.
    void construct_group_dof(Space fs, const Key group_key={0,0}, int p_order=1);


};

}