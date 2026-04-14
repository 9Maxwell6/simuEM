#pragma once

#include "utils/logger.h"
#include "utils/util_la.h"

#include "math/fem/block.h"
#include "math/data_format.h"




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
    friend class FEM_System;
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

    std::string print_block_rack() const;

    void delete_data();

    void build_linear_system();



    //TODO
    //void initial_guess();



    



};


}