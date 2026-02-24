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
    size_t dof_offset_;

    std::vector<size_t> dof_space_offset_;   // used for block matrix assemble.

    std::vector<size_t> elements_dof_lookup_list;  // [dofs of element 1, dofs of element 2, ...]

    std::vector<FEM_Space * > global_space;
    std::unordered_map<Key, std::vector<FEM_Space * >, Key_Hash> group_space;


    // multiple fe_space will result in block matrices configuration; 
    // hence, we adopt a convention for block positions:
    //  
    //   dof_space_offset_ = [ 0,               ...,  global_n,        group_1,        ...,  group_m         ]
    //                         |                 |    |                |                |    |
    //                         |                 |    |                |                |    |
    //                         ↓                 |    |                |                |    |
    //                       [ [global_space_1]  ↓.   |      .         |     .          |.   |      .        ] 
    //                       [        .          ...  ↓      .         |     . {coupling|between spaces}     ]
    //                       [        .           .   [global_space_n] ↓     .          |.   |     .         ]
    //                       [        .           .          .         [group_space_1]  ↓.   |     .         ]
    //                       [   {coupling between spaces}   .               .          ...  ↓     .         ]
    //                       [        .           .          .               .           .   [group_space_m] ]
    //
    //
    //
    // for each block space, dof is ordered from:
    //      node -> edge -> face -> volume
    // i.e.:
    //
    //
    //      [ [nodes]    .       .        .     ]
    //      [    .    [edges]    .        .     ]
    //      [    .       .    [faces]     .     ]
    //      [    .       .       .    [volumes] ]
    //
    //




    //std::vector<size_t> elements_dof_lookup_list;

    //void assign_dof(Element * e);

    bool initialize_space_dof();

    size_t assign_element_dof(FEM_Space& fe_space, Element& e);


    // assume conforming FEM mesh
    std::array<size_t, 4> count_node_edge_face_volume_dof(const std::vector<Element*>& elements);

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