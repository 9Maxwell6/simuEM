#pragma once

#include "math/fem/block.h"
#include "math/data_format.h"

#include "utils/logger.h"
#include "utils/util_la.h"

// for test only
//#include "utils/util_petsc.h"



namespace simu {

/**
 * 
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 */
class Block_Rack
{
    friend class  FEM_System;
    friend struct Dirichlet_BC;
private:
    
    size_t n_row_; 
    size_t n_col_; 

    std::vector<Block*>  rack_;

    std::vector<size_t> unit_row_length_;
    std::vector<size_t> unit_col_length_;

    G_Matrix lhs_;
    G_Vector rhs_;
    G_Vector   x_;    // solution

    
public:
    Block_Rack(){};
    Block_Rack(size_t n_row, size_t n_col);

    void set_grid(size_t n_row, size_t n_col);

    bool insert_block(Block& block, size_t row, size_t col);
    bool compute_block_offset();

    bool is_block_ready(const Block* block) const;

    void build_linear_system();


    void delete_data();
    

    std::string print_block_rack() const;


    size_t get_n_row() const { return n_row_; }
    size_t get_n_col() const { return n_col_; }

    const G_Matrix get_lhs() const { return lhs_; }
    const G_Vector get_x()   const { return x_;   }
    const G_Vector get_rhs() const { return rhs_; }

    const std::vector<Block*>& get_rack() const { return rack_; }
    //TODO:
    //void initial_guess();

    //TODO: Dirichlet BC

    

};


}