#pragma once

#include "math/fem/assemble_data.h"
#include "math/fem/element_data.h"
#include "math/fem/shape.h"
#include "math/operator/operator_collection.h"




namespace simu {


template<typename Op>
bool assemble_vec(const Assemble_Data& data, Op&& user_operation)
{
    
    if(data.block_vector == nullptr)
    {
        Logger::error("assemble_vec: global vector not initialized, call FEM_System::assemble_vec_data first.");
        return 0;
    }

    auto assemble_local = [&](auto& e_data, auto& local_vec, const auto& user_operation)
    {        
        local_vec.setZero();

        user_operation(e_data, local_vec);

        Operator::dof_transformation_vec(e_data, local_vec);

        Operator::add_to_global_vec(data, local_vec);

    };


    auto assemble_global = [&](auto& e_data, const auto& user_operation) 
    {
        e_data.mesh = data.mesh;

        e_data.space_1 =  data.space_1;

        int vdim_1 = data.space_1->get_vdim();

        e_data.check = &data.check;

        data.row_dof_offset = 0;
        data.col_dof_offset = 0;

        for(const Element* e : *data.elements)
        {
            
            //int order = e->get_geometry_order() + space_1->get_basis_order() + space_2->get_basis_order();
            Basis_Shape b_shape = to_basis_shape(e->get_geometry());

            e_data.e = e;
            e_data.b_shape = b_shape;

            e_data.shape_space_1 = data.space_1->get_basis_space(b_shape);

            size_t row_size = e_data.shape_space_1->get_n_dof();

            size_t e_row = row_size * vdim_1;

            e_data.rows = e_row;

            e_data.reset_flag();


            // optimized with fixed size matrix for certain cases.
            if      (e_row == 2 ) { Vector<2> vec;                     assemble_local(e_data, vec, user_operation); }
            else if (e_row == 3 ) { Vector<3> vec;                     assemble_local(e_data, vec, user_operation); }
            else if (e_row == 4 ) { Vector<4> vec;                     assemble_local(e_data, vec, user_operation); }
            else if (e_row == 6 ) { Vector<6> vec;                     assemble_local(e_data, vec, user_operation); }
            else                  { VectorXd  vec(e_row);              assemble_local(e_data, vec, user_operation); }

            data.row_dof_offset += row_size;
        }
    };

    int phy_dim = data.mesh_dim;
    int ref_dim = data.element_dim;

    if      (phy_dim == 3 && ref_dim == 3) { static Element_Data<3,3> e_data; assemble_global(e_data, user_operation); }
    else if (phy_dim == 3 && ref_dim == 2) { static Element_Data<3,2> e_data; assemble_global(e_data, user_operation); }
    else if (phy_dim == 3 && ref_dim == 1) { static Element_Data<3,1> e_data; assemble_global(e_data, user_operation); }
    else if (phy_dim == 2 && ref_dim == 2) { static Element_Data<2,2> e_data; assemble_global(e_data, user_operation); }
    else if (phy_dim == 2 && ref_dim == 1) { static Element_Data<2,1> e_data; assemble_global(e_data, user_operation); }
    else if (phy_dim == 1 && ref_dim == 1) { static Element_Data<1,1> e_data; assemble_global(e_data, user_operation); }


    la_kernel::finalize_vec(data.block_vector);

    return true;
}



}