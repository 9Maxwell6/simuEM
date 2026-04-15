#include "math/fem/block.h"


using namespace simu;


void Block::block_transpose(const Block& block)
{
    row_size = block.col_size;
    col_size = block.row_size;

    if(row_size != block.col_size) 
    {
        Logger::warning("Block::block_transpose: dangeous - "
                        "The block's row_size [id="+std::to_string(id)+", row_size="+std::to_string(row_size)+"] doesn't match "
                        "the target block's col_size [id="+std::to_string(block.id)+", col_size="+std::to_string(block.col_size)+"], "
                        "so they've been forced to match. Please avoid using assemble_mat on this block in the future.");
        row_size = block.col_size;
    }

    if(col_size != block.row_size) 
    {
        Logger::warning("Block::block_transpose: dangeous - "
                        "The block's col_size [id="+std::to_string(id)+", row_size="+std::to_string(col_size)+"] doesn't match "
                        "the target block's row_size [id="+std::to_string(block.id)+", col_size="+std::to_string(block.row_size)+"], "
                        "so they've been forced to match. Please avoid using assemble_mat on this block in the future.");
        col_size = block.row_size;
    }

    if(is_base_block != block.is_base_block) 
    {
        Logger::warning("Block::block_transpose: unusual case - is_base_block flag not match: "
                        "[id="+std::to_string(id)+", is_base_block="+std::to_string(is_base_block)+"] and "
                        "[id="+std::to_string(block.id)+", is_base_block="+std::to_string(block.is_base_block)+"]");
    }

    if(block.is_base_block) Logger::warning("Block::block_transpose: unusual case - trying to transpose target block at diagonal entry.");


    // destory original matrix
    la_kernel::destroy_mat(mat);
    la_kernel::create_transpose(block.mat, mat);

}