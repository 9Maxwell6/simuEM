#pragma once

#include "math/fem/space_collection.h"
#include "entity/mesh/e_collection.h"
#include "math/fem/dof_handler.h"
#include "entity/mesh/mesh.h"


#include "utils/util_hash.h"



namespace simu {


struct Block 
{
    size_t id; 
    size_t row_offset;
    size_t col_offset;
    size_t row_size;
    size_t col_size;

    bool operator==(const Block& other) const 
    {
        return id == other.id;
    }

    struct Hash 
    {
        size_t operator()(const Block& b) const 
        {
            return std::hash<size_t>{}(b.id);
        }
    };
};


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
    DoF_Handler dof_handler_;

    std::vector<size_t> dof_list_;     // entry -> index in global dof
    size_t dof_offset_;                // for creating next dof

    std::vector<size_t> dof_space_list_;   // used for block matrix assemble. (starting indices in the global dof list)
    size_t dof_space_offset_;              // for creating next dof for the block matrix

    std::vector<size_t> elements_dof_lookup_list;  // [dofs of element 1, dofs of element 2, ...]


    std::vector<FEM_Space * > global_space_;


    std::unordered_map<Key, std::vector<FEM_Space * >, Key::Hash> group_space_;

    
    size_t block_id_;
    std::unordered_map<Block, FEM_Space *,         Block::Hash> fe_block_space_;
    std::unordered_map<Block, Key,                 Block::Hash> fe_block_key_;
    std::unordered_map<Block, std::vector<size_t>, Block::Hash> fe_block_dof_;

    std::unordered_map<Block, std::pair<Block, Block>, Block::Hash> coupled_block_;




    


    // coupling between space need two block:   and dof list from each block



    //std::vector<size_t> elements_dof_lookup_list;

    //void assign_dof(Element * e);

    bool generate_block_dof(Block& block);

    bool generate_coupling_block_dof(Block& block_1, Block& block_2);





    static Basis_Shape to_basis_shape(Geometry t);

    static Geometry to_element_geometry(Basis_Shape t);


public:
    FEM_System(Mesh& mesh);


    // assign functional space to specific group of elements, 
    // if using default key, assign space to global domain.
    Block register_FE_space(FEM_Space& fe_space, const Key group_key={0,0});

    Block register_FE_space_coupling(const Block& block_1, const Block& block_2, const Key group_key={0,0});


    
    const FEM_Space* get_block_space(const Block& block) const;
    const Key get_block_group_key(const Block& block) const;
    const std::vector<size_t>& get_block_dof(const Block& block) const;

    const std::pair<Block, Block>& get_coupled_block(const Block& block) const;


};

}