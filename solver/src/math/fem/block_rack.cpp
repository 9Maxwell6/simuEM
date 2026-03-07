#include "math/fem/block_rack.h"

using namespace simu;

Block_Rack::Block_Rack(size_t n_row, size_t n_col) : n_row_(n_row), n_col_(n_col)
{
    rack_.resize(n_row * n_col);
    unit_row_length_.resize(n_row);
    unit_col_length_.resize(n_col);
}

void Block_Rack::set_grid(size_t n_row, size_t n_col)
{
    n_row_ = n_row;
    n_col_ = n_col;
    rack_.resize(n_row * n_col);
    unit_row_length_.resize(n_row);
    unit_col_length_.resize(n_col);
}

bool Block_Rack::insert_block(Block& block, size_t row, size_t col)
{
    if(row>=n_row_ || col>=n_col_) {Logger::error("Block_Rack::insert_block - invalid rack slot."); return false; }
    if(rack_[row*n_col_ +col]!= nullptr) Logger::warning("Block_Rack::insert_block - slot already occupied, replaing with new block."); 

    if(unit_row_length_[row]!=0 && block.row_size!=unit_row_length_[row]) {Logger::error("Block_Rack::insert_block - block row size mismatch."); return false; }
    if(unit_col_length_[col]!=0 && block.col_size!=unit_row_length_[row]) {Logger::error("Block_Rack::insert_block - block row size mismatch."); return false; }

    unit_row_length_[row] = block.row_size;
    unit_col_length_[col] = block.col_size;

    rack_[row*n_col_ +col] = &block;

    return true;
}

bool Block_Rack::compute_block_offset()
{
    std::vector<size_t> row_offsets(n_row_);
    std::vector<size_t> col_offsets(n_col_);

    row_offsets[0] = 0;
    for (int i = 1; i < n_row_; ++i)
        row_offsets[i] = row_offsets[i - 1] + unit_row_length_[i - 1];

    col_offsets[0] = 0;
    for (int j = 1; j < n_col_; ++j)
        col_offsets[j] = col_offsets[j - 1] + unit_col_length_[j - 1];

    for (int i = 0; i < n_row_; ++i) {
        for (int j = 0; j < n_col_; ++j) {
            Block * b = rack_[i*n_col_+j];
            if(b==nullptr) {Logger::error("Block_Rack::compute_block_offset - block rack has empty slot."); return false; }
            b->row_offset = row_offsets[i];
            b->col_offset = col_offsets[j];
        }
    }

    return true;
}