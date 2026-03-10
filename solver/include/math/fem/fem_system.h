#pragma once

#include "math/fem/space_collection.h"
#include "entity/mesh/e_collection.h"
#include "math/fem/block_rack.h"
#include "entity/mesh/mesh.h"


#include "utils/util_hash.h"

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
 *     node -> edge -> face -> volume
 * i.e.:
 *
 *      [ [nodes]    .       .        .     ]
 *      [    .    [edges]    .        .     ]
 *      [    .       .    [faces]     .     ]
 *      [    .       .       .    [volumes] ]
 *
 * the mapping between 
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


    std::vector<FEM_Space * > global_space_;


    std::unordered_map<Key, std::vector<FEM_Space * >, Key::Hash> group_space_;

    
    size_t block_id_;
    std::unordered_map<Block, Key,                 Block::Hash> fe_block_key_;

    // for basic block
    std::unordered_map<Block, FEM_Space *,         Block::Hash> fe_block_space_;
    std::unordered_map<Block, const std::vector<size_t> *, Block::Hash> fe_block_dof_;
    //temporary, for constructing dof, should be cleared after dof table is constructed.
    std::unordered_map<Block, util::Block_Hash,    Block::Hash> fe_block_hash_;

    // for coupling block
    std::unordered_map<Block, std::array<Block, 2>,                       Block::Hash> coupled_block_;
    std::unordered_map<Block, std::array<FEM_Space *, 2>,                 Block::Hash> coupled_block_space_;
    std::unordered_map<Block, std::array<const std::vector<size_t> *, 2>, Block::Hash> coupled_block_dof_;

    // store actual coupling block dof data, this is to avoid copy when transpose of block is applied.
    std::unordered_map<Block, std::vector<size_t>,   Block::Hash> fe_block_dof_data_;
    std::unordered_map<Block, std::array<std::vector<size_t>, 2>,   Block::Hash> coupled_block_dof_data_;




    


    // coupling between space need two block:   and dof list from each block



    //std::vector<size_t> elements_dof_lookup_list;

    //void assign_dof(Element * e);

    const util::Block_Hash_D create_volume_dof_hash(const Block& block) const;

    bool generate_block_dof(Block& block);

    bool generate_coupling_block_dof(Block& block);



    static Basis_Shape to_basis_shape(Geometry t);

    static Geometry to_element_geometry(Basis_Shape t);


    template <typename Get_dof>
    bool node_edge_face_dof_handler(Get_dof&& dof_handler, Basis_Shape shape, const size_t* node_idx, int node_size, 
                                                                    int n_dof_per_node, size_t node_dof_offset,
                                                                    int n_dof_per_edge, size_t edge_dof_offset, 
                                                                    int n_dof_per_face, size_t face_dof_offset);


public:
    FEM_System(Mesh& mesh);


    // assign functional space to specific group of elements, 
    // if using default key, assign space to global domain.
    Block register_FE_space(FEM_Space& fe_space, const Key group_key={0,0});

    Block register_FE_space_coupling(const Block& block_1, const Block& block_2, const Key group_key={0,0});


    
    const FEM_Space* get_block_space(const Block& block) const;
    const Key get_block_group_key(const Block& block) const;
    const std::vector<size_t> * get_block_dof(const Block& block) const;
    const util::Block_Hash& get_block_hash(const Block& block) const;

    const std::array<Block, 2>& get_coupled_block(const Block& block) const;
    const std::array<FEM_Space *, 2>& get_coupled_block_space(const Block& block) const;
    const std::array<const std::vector<size_t> * ,2>& get_coupled_block_dof(const Block& block) const; 

    // use only after every dof table of blocks are initialized
    void delete_block_hash() { fe_block_hash_.clear(); }

    // TODO:
    Block transpose_block(const Block& block);

    Block_Rack initialize_block_rack(size_t n_row, size_t n_col);

};

}