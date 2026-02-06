#pragma once
#include "element.h"

#include <vector>

class Mesh
{
protected:
    int dim_; // space dimension
    size_t n_node    , n_edge    , n_face    ;  // actual number of node/edge/face in mesh
    //size_t N_dof_node, N_dof_edge, N_dof_face; 
    size_t n_element;

    size_t n_exterior_boundary_node; // true boundary of simulation domain
    size_t n_interior_boundary_node;

    std::vector<Element> elements;

    



};