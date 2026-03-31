#pragma once

#include "math/fem/assemble_data.h"
#include "math/fem/integrator.h"
#include "math/fem/operation.h"
#include "math/fem/fem_util.h"



namespace simu {



template<typename Op>
bool assemble_block(const Assemble_Data& data, Op&& user_operation)
{
    

    auto assemble_local = [&](auto& e_data, auto& local_mat, const auto& user_operation)
    {        
        local_mat.setZero();

        user_operation(e_data, local_mat);

        Operation::dof_transformation(e_data, local_mat);


    };


    auto assemble_global = [&](auto& e_data, const auto& user_operation) 
    {
        size_t col_dof_counter = 0;
        size_t row_dof_counter = 0;

        e_data.mesh = data.mesh;

        e_data.space_1 =  data.space_1->get_function_space();
        e_data.space_2 =  data.space_2->get_function_space();

        e_data.integrator_check = &data.integrator_check_flags;

        for(const Element* e : *data.elements)
        {
            
            //int order = e->get_geometry_order() + space_1->get_basis_order() + space_2->get_basis_order();
            Basis_Shape b_shape = to_basis_shape(e->get_geometry());

            e_data.e = e;
            e_data.b_shape = b_shape;
            e_data.i_r_list = &data.integration_rule[b_shape];

            e_data.shape_space_1 = data.space_1->get_basis_space(b_shape);
            e_data.shape_space_2 = data.space_2->get_basis_space(b_shape);

            size_t row_size = e_data.shape_space_1->get_n_dof();
            size_t col_size = e_data.shape_space_2->get_n_dof();

            e_data.rows = row_size;
            e_data.cols = col_size;

            e_data.reset_flag();


            // optimized with fixed size matrix for certain cases.
            if      (row_size == 3 && col_size == 3) { Matrix<3,3> mat;                     assemble_local(e_data, mat, user_operation); }
            else if (row_size == 4 && col_size == 4) { Matrix<4,4> mat;                     assemble_local(e_data, mat, user_operation); }
            else if (row_size == 4 && col_size == 6) { Matrix<4,6> mat;                     assemble_local(e_data, mat, user_operation); }
            else if (row_size == 6 && col_size == 4) { Matrix<6,4> mat;                     assemble_local(e_data, mat, user_operation); }
            else if (row_size == 6 && col_size == 6) { Matrix<6,6> mat;                     assemble_local(e_data, mat, user_operation); }
            else                                     { MatrixXd    mat(row_size, col_size); assemble_local(e_data, mat, user_operation); }

            row_dof_counter += row_size;
            col_dof_counter += col_size;
            
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



    return true;
}





}