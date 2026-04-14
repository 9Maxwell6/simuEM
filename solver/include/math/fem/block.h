#pragma once

#include "utils/logger.h"
#include "math/data_format.h"



namespace simu {

struct Block 
{
    size_t id; 
    size_t row_offset;
    size_t col_offset;
    size_t row_size;
    size_t col_size;

    bool is_base_block;

    G_Matrix mat;
    G_Vector vec;  // only initialized if is_base_block == true; 


    bool operator==(const Block& other) const { return id == other.id; }

    struct Hash { size_t operator()(const Block& b) const { return std::hash<size_t>{}(b.id);} };
};



}