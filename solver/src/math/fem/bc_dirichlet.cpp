#include "math/fem/bc_dirichlet.h"

using namespace simu;


bool Dirichlet_BC::apply_to_system(Block_Rack& br)
{
    bool match_block = false;
    size_t row_offset = 0;
    for(Block* br_block : br.rack_) 
    {
        if(block == br_block){
            match_block = true;
            row_offset = br_block->row_offset;
        }
    }

    if(!match_block) {Logger::error("Dirichlet_BC::apply_to_system: block rack does not contain the block attached to this Dirichlet BC."); return false;}
    
    switch (bc_type) 
    {
        case Dirichlet_Type::HOMOGENEOUS: { la_kernel::set_value_vec(bc_dofs, std::vector<scalar_t>(bc_dofs.size(), 0.), br.x_); break; }
        case Dirichlet_Type::CONSTANT:    { la_kernel::set_value_vec(bc_dofs, std::vector<scalar_t>(bc_dofs.size(), bc_values[0]), br.x_); break;}
        case Dirichlet_Type::FIELD:       { la_kernel::set_value_vec(bc_dofs, bc_values, br.x_); break;}
    }

    std::vector<dof_idx> offset_dofs(bc_dofs.size());
    for (dof_idx i = 0; i < bc_dofs.size(); ++i) offset_dofs[i] = bc_dofs[i] + row_offset;


    // TODO: consider vdim and layout for vector H1
    la_kernel::zero_row_col_mat(offset_dofs, 1., br.lhs_, br.x_, br.rhs_);


    return true;
}


bool Dirichlet_BC::apply_to_system(G_Matrix lhs, G_Vector rhs, G_Vector x)
{
    switch (bc_type) 
    {
        case Dirichlet_Type::HOMOGENEOUS: { la_kernel::set_value_vec(bc_dofs, std::vector<scalar_t>(bc_dofs.size(), 0.), x); break; }
        case Dirichlet_Type::CONSTANT:    { la_kernel::set_value_vec(bc_dofs, std::vector<scalar_t>(bc_dofs.size(), bc_values[0]), x); break;}
        case Dirichlet_Type::FIELD:       { la_kernel::set_value_vec(bc_dofs, bc_values, x); break;}
    }

    // TODO: consider vdim and layout for vector H1
    la_kernel::zero_row_col_mat(bc_dofs, 1., lhs, x, rhs);
    return true;
}