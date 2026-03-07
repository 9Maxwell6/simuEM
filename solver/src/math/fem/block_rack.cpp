#include "math/fem/block_rack.h"

using namespace simu;

Block_Rack::Block_Rack(size_t row_block_size, size_t col_block_size) : row_block_size_(row_block_size), col_block_size_(col_block_size)
{
    rack_.resize(row_block_size*col_block_size);
    unit_row_length_.resize(row_block_size);
    unit_col_length_.resize(col_block_size);
}

void Block_Rack::set_grid(size_t row_block_size, size_t col_block_size)
{
    rack_.resize(row_block_size*col_block_size);
    unit_row_length_.resize(row_block_size);
    unit_col_length_.resize(col_block_size);
}

bool Block_Rack::insert_block(const Block& block, size_t row, size_t col)
{
    if(rack_[row*])
}