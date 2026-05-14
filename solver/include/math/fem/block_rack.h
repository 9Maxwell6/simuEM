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

    std::vector<size_t> unit_row_length_;  // global rows in each block-row
    std::vector<size_t> unit_col_length_;  // global cols in each block-col

    // In a serial run (MPI comm size == 1), local_row/col_size_ == unit_row/col_length_
    std::vector<size_d> local_row_size_;   // local rows in each block-row on this rank
    std::vector<size_d> local_col_size_;   // local cols in each block-col on this rank

    G_Matrix lhs_ = nullptr;
    G_Vector rhs_ = nullptr;
    G_Vector   x_ = nullptr;    // solution


    // convert global big matrix/vector to block matrices and vectors, used for block operations.
    std::vector<G_Matrix>  block_lhs_;
    std::vector<G_Vector>  block_rhs_;
    std::vector<G_Vector>  block_x_;

    

    
public:
    Block_Rack(){};
    Block_Rack(size_t n_row, size_t n_col);
    ~Block_Rack();

    void set_grid(size_t n_row, size_t n_col);

    bool insert_block(Block& block, size_t row, size_t col);
    bool compute_block_offset();

    bool is_block_ready(const Block* block) const;

    void build_linear_system();

    void extract_block_system();


    void delete_data();
    

    std::string print_block_rack() const;


    size_t get_n_row() const { return n_row_; }
    size_t get_n_col() const { return n_col_; }

    const G_Matrix get_lhs() const { return lhs_; }
    const G_Vector get_x()   const { return x_;   }
    const G_Vector get_rhs() const { return rhs_; }
    
    // must call extract_block_system() before get_mat or get_vec
    const G_Matrix get_block_lhs(size_t row_idx, size_t col_idx) const { return block_lhs_[row_idx*n_col_+col_idx]; }
    const G_Vector get_block_rhs(size_t row_idx                ) const { return block_rhs_[row_idx];                }
    const G_Vector get_block_x  (size_t row_idx                ) const { return block_x_[row_idx];                  }

    void assemble_block_lhs();
    void assemble_block_rhs();
    void assemble_block_x();


    const std::vector<Block*>& get_rack() const { return rack_; }
    //TODO:
    //void initial_guess();

    //TODO: Dirichlet BC

    void finalize();

    

};


}