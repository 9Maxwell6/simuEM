#pragma once

#include "utils/logger.h"


#include <cstddef>
#include <functional>


namespace simu {

struct Block 
{
    size_t id; 
    size_t row_offset;
    size_t col_offset;
    size_t row_size;
    size_t col_size;

    bool is_base_block;


    bool operator==(const Block& other) const { return id == other.id; }

    struct Hash { size_t operator()(const Block& b) const { return std::hash<size_t>{}(b.id);} };
};



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

    std::vector<Block*> rack_;

    std::vector<size_t> unit_row_length_;
    std::vector<size_t> unit_col_length_;

    
public:
    Block_Rack(){};
    Block_Rack(size_t n_row, size_t n_col);

    void set_grid(size_t n_row, size_t n_col);

    bool insert_block(Block& block, size_t row, size_t col);

    bool compute_block_offset();

    std::string print_block_rack() const;


    



};


}