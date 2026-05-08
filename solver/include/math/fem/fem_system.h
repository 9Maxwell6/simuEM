#pragma once

#include "config.h"

#include "math/fem/space_collection.h"
#include "world/mesh/e_collection.h"
#include "math/fem/block_rack.h"
#include "world/mesh/mesh.h"
#include "math/data_format.h"
#include "math/fem/assemble_data.h"
#include "math/fem/shape.h"

#include "math/fem/bc_dirichlet.h"


//#include "world/mesh/e__transformation.h"



#include "utils/util_hash.h"
#include "utils/util_petsc.h"
#include "utils/util_la.h"


#include <utility>



namespace simu {



/**
 * multiple fe_space will result in block matrices configuration; 
 * hence, we adopt a convention for block positions:
 *  
 *     dof_space_list_ = [ 0,               ...,  global_n,        group_1,        ...,  group_m         ]
 *                         |                 |    |                |                |    |
 *                         |                 |    |                |                |    |
 *                         ↓                 |    |                |                |    |
 *                       [ [global_space_1]  ↓.   |      .         |     .          |.   |      .        ] 
 *                       [        .          ...  ↓      .         |     . {coupling|between spaces}     ]
 *                       [        .           .   [global_space_n] ↓     .          |.   |     .         ]
 *                       [        .           .          .         [group_space_1]  ↓.   |     .         ]
 *                       [   {coupling between spaces}   .               .          ...  ↓     .         ]
 *                       [        .           .          .               .           .   [group_space_m] ]
 *
 *
 *
 * for each block space, dof is ordered from:
 *     node -> edge -> face -> cell
 * i.e.:
 *
 *      [ [nodes]    .       .        .     ]
 *      [    .    [edges]    .        .     ]
 *      [    .       .    [faces]     .     ]
 *      [    .       .       .    [cells] ]
 *
 * the mapping between 
 * 
 * TODO: separate dof logic from FEM_System, move all dof logic into math/dof folder.
 */
class FEM_System
{

private:
    int dim_;

    Mesh& mesh_;
    Block_Rack block_rack_;

    std::vector<size_t> dof_list_;     // entry -> index in global dof
    size_t dof_offset_;                // for creating next dof

    std::vector<size_t> dof_space_list_;   // used for block matrix assemble. (starting indices in the global dof list)
    size_t dof_space_offset_;              // for creating next dof for the block matrix

    std::vector<size_t> elements_dof_lookup_list;  // [dofs of element 1, dofs of element 2, ...]


    //std::vector<FEM_Space * > global_space_;


    //std::unordered_map<Key, std::vector<FEM_Space * >, Key::Hash> group_space_;

    
    size_t block_id_;
    std::unordered_map<Block, Key,                          Block::Hash> fe_block_key_;

    // for basic block
    std::unordered_map<Block, FEM_Space *,                  Block::Hash> fe_block_space_;
    std::unordered_map<Block, const std::vector<dof_idx> *, Block::Hash> fe_block_dof_;
    std::unordered_map<Block, util::Block_Hash,             Block::Hash> fe_block_hash_;

    // for coupling block
    std::unordered_map<Block, std::array<const Block*, 2>,                  Block::Hash> coupled_block_;
    std::unordered_map<Block, std::array<const FEM_Space*, 2>,              Block::Hash> coupled_block_space_;
    std::unordered_map<Block, std::array<const std::vector<dof_idx>*, 2>,   Block::Hash> coupled_block_dof_;
    std::unordered_map<Block, std::array<util::Block_Hash, 2>,              Block::Hash> coupled_block_hash_;

    // store actual coupling block dof data, this is to avoid copy when transpose of block is applied.
    std::unordered_map<Block, std::vector<dof_idx>,   Block::Hash> fe_block_dof_data_;
    std::unordered_map<Block, std::array<std::vector<dof_idx>, 2>,   Block::Hash> coupled_block_dof_data_;

    

    // coupling between space need two block:   and dof list from each block



    //std::vector<size_t> elements_dof_lookup_list;

    //void assign_dof(Element * e);

    const util::Block_Hash_D create_cell_dof_hash(const Block& block) const;

    bool generate_block_dof(Block& block, FEM_Space& fe_space, int idx);

    bool generate_coupling_block_dof(Block& block);

    // compute number of non-zero entry per row, used for block matrix pre-allocation.                                                                
    std::vector<size_d> compute_nnz_per_row(
        const FEM_Space* space_1, const std::vector<dof_idx>* block_row_dof, size_d block_row_size, 
        const FEM_Space* space_2, const std::vector<dof_idx>* block_col_dof, size_d block_col_size,
        const std::vector<Element*>* elements) const;


    template <typename Get_dof>
    bool node_edge_face_dof_handler(Get_dof&& dof_handler, Basis_Shape shape, const size_t* node_idx, int node_size, 
                                                                    int n_dof_per_node, size_t node_dof_offset,
                                                                    int n_dof_per_edge, size_t edge_dof_offset, 
                                                                    int n_dof_per_face, size_t face_dof_offset);
    

    /*
    // not used
    // used for transformation from reference element to actual element in mesh.
    Transform_general_2D          transform_element_2D;
    Transform_general_3D          transform_element_3D;

    // special optimization is implemented:
    Transform_Triangle_o1_2D      transform_triangle_o1_2D;
    Transform_Triangle_o1_3D      transform_triangle_o1_3D;
    Transform_Tetrahedron_o1_3D   transform_tetrahedron_o1_3D;
    */


public:
    FEM_System(Mesh& mesh);


    // assign functional space to specific group of elements, 
    // if using default key, assign space to global domain.
    Block register_FE_space(FEM_Space& fe_space, const Key group_key={0,0}, const Block* block=nullptr);

    Block register_dual_FE_space(FEM_Space& fe_space_1, FEM_Space& fe_space_2, const Key group_key={0,0}, const Block* block_1=nullptr, const Block* block_2=nullptr);

    Block register_FE_space_coupling(const Block& block_1, const Block& block_2, const Key group_key={0,0});

    Dirichlet_BC register_Dirichlet_BC(const Block& block, const Key& group_key, Dirichlet_Type bc_type);

    // important !!!!
    // TODO: when assign dof, for each block, record dof of boundary, separate by group key.
    // [dof of element 1, dof of element 2]

    const FEM_Space* get_block_space(const Block& block) const;
    const Key get_group_key(const Block& block) const;
    const std::vector<dof_idx> * get_block_dof(const Block& block) const;
    const util::Block_Hash& get_block_hash(const Block& block) const;

    const std::array<const Block*,                 2>& get_coupled_block(const Block& block) const;
    const std::array<const FEM_Space *,            2>& get_coupled_block_space(const Block& block) const;
    const std::array<const std::vector<dof_idx> * ,2>& get_coupled_block_dof(const Block& block) const; 

    const FEM_Space*            get_block_space(const Block& block, int idx) const;
    const std::vector<dof_idx>* get_block_dof(const Block& block, int idx) const; 

    const FEM_Space* get_block_row_space(const Block& block) const;
    const FEM_Space* get_block_col_space(const Block& block) const;
    const std::vector<dof_idx>* get_block_row_dof(const Block& block) const;
    const std::vector<dof_idx>* get_block_col_dof(const Block& block) const;

    // use only after every dof table of blocks are initialized
    void delete_block_hash() { fe_block_hash_.clear(); }

    const Mesh& get_mesh() const { return mesh_; }


    // TODO:
    Block transpose_block(const Block& block);

    Block_Rack initialize_block_rack(size_t n_row, size_t n_col);

    
    Assemble_Data assemble_mat_data(Block& block);
    
    Assemble_Data assemble_vec_data(Block& block);  // block must be base_block.



};

}