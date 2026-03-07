#pragma once
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
    
    size_t row_block_size_; 
    size_t col_block_size_; 

    std::vector<Block*> rack_;

    std::vector<size_t> unit_row_length_;
    std::vector<size_t> unit_col_length_;
    
public:
    Block_Rack(){};
    Block_Rack(size_t row_block_size, size_t col_block_size);

    void set_grid(size_t row_block_size, size_t col_block_size);

    bool insert_block(const Block& block, size_t row, size_t col);





};


}