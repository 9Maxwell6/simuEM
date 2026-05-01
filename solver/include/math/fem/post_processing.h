#pragma once

#include "math/fem/fem_system.h"
#include "math/fem/block_rack.h"
#include "math/fem/shape.h"
#include "math/data_format.h"

#include "utils/util_hash.h"
#include "utils/util_la.h"

#include <unordered_set>
#include <tuple>

namespace simu {


struct Element_DoF {
    const Block* block;

    size_t dof_offset;
    size_t dof_size;
};


template<typename Op>
scalar_t integrate_element(const Block_Rack& br_system, const FEM_System& fe_system, Op&& user_operation)
{

    auto global_integral = [&](auto& e_data, const Block_Rack& br_system, const FEM_System& fe_system, const auto& user_operation)
    {        
       
        const Mesh& mesh = fe_system.get_mesh();

        std::vector<scalar_t> x;
        la_kernel::extract_vec(br_system.get_x(), x);

        scalar_t result = 0;

        

        const std::vector<Block*>& rack = br_system.get_rack();

        size_t n_row = br_system.get_n_row();
        size_t n_col = br_system.get_n_col();


        // map between element and the column index of the base block containing it. 
        // e.g. element is contained in base block 0, 1, 4.  {e*, {0,1,4}}
        std::unordered_set<Element*> e_record;


        for(size_t i=0; i<n_row; ++i)
        {
            const FEM_Space* space_row;

            // map between element and column index of the block containing the element.
            std::unordered_map<Element*, std::vector<Element_DoF>> element_map;

            for(size_t j=0; j<n_col; ++j)
            {
                const Block* block = rack[i*n_row+j];

                const Key& group_key = fe_system.get_block_group_key(*block);
                const std::vector<Element*>& elements = mesh.get_element_group(group_key);
                
                size_t col_dof_offset = 0;
                for(Element* e : elements)
                {
                    Basis_Shape b_shape = to_basis_shape(e->get_geometry());
                    const FEM_Space* shape_space = fe_system.get_block_col_space(*block)->get_basis_space(b_shape);

                    size_t col_size = shape_space->get_n_dof();

                    const std::vector<dof_idx>* col_dof_list = fe_system.get_block_col_dof(*block);

                    if(e_record.count(e) == 0){
                        // new element
                        std::vector<Element_DoF>& e_dof_list = element_map[e];
                        Element_DoF e_dof{.block = block, .dof_offset = col_dof_offset, .dof_size = col_size};
                        //e_dof.block = block;
                        //e_dof.dof_offset = col_dof_offset;
                        //e_dof.dof_size = col_size;
                        e_dof_list.push_back(e_dof);
                    }
                    col_dof_offset += col_size;
                }

            }

            for(size_t j=0; j<n_col; ++j)
                for(Element* e : mesh.get_element_group(fe_system.get_block_group_key(*rack[i*n_row+j]))) e_record.insert(e);
            

            int ccc = 0;
            for (const auto& [e, e_dof_list] : element_map) {

                e_data.mesh = &mesh;
                e_data.e = e;

                e_data.b_shape = to_basis_shape(e->get_geometry());

                e_data.reset_flag();


                std::vector<const FEM_Space*> space_list;
                std::vector<std::vector<scalar_t>> dof_value_list;

                for(const Element_DoF& e_dof : e_dof_list) if(e_dof.block)
                {

                    //const FEM_Space* row_space = fe_system.get_block_row_space(*block);
                    //const FEM_Space* col_space = fe_system.get_block_col_space(*block);
                    //const std::vector<dof_idx>* row_dof_list = fe_system.get_block_row_dof(*block);
                    const std::vector<dof_idx>* col_dof_list = fe_system.get_block_col_dof(*e_dof.block);

                    space_list.push_back(fe_system.get_block_col_space(*e_dof.block)->get_basis_space(e_data.b_shape));

                    std::vector<scalar_t> block_dof_value;
                    for(size_t idx=0; idx<e_dof.dof_size; ++idx)
                    {
                        block_dof_value.push_back(x[e_dof.block->col_offset + (*col_dof_list)[e_dof.dof_offset + idx]]);
                    }
                    dof_value_list.push_back(std::move(block_dof_value));
                }

                e_data.space_list = &space_list;
                e_data.dof_value_list = &dof_value_list;
                
                user_operation(e_data, result);

            }
        }
        return result;

    };



    
    int phy_dim = fe_system.get_mesh().get_mesh_dimension();


    // TODO: extend to 1D and 2D
    if      (phy_dim == 3) { static Element_Data<3,3> e_data; return global_integral(e_data, br_system, fe_system, user_operation); }
    //else if (phy_dim == 2) { static Element_Data<2,2> e_data; return global_integral(e_data, br_system, fe_system, user_operation); }
    //else if (phy_dim == 1) { static Element_Data<1,1> e_data; return global_integral(e_data, br_system, fe_system, user_operation); }

    return 0.;

}


}