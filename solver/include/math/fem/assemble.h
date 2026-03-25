#pragma once

#include "math/fem/assemble_data.h"
#include "math/fem/integrator.h"
#include "math/fem/fem_util.h"



namespace simu {



template<typename Op>
bool assemble_block(const Assemble_Data& data, Op&& user_operation)
{
    

    auto assemble_local = [&](const Element* e, auto& local_mat, const auto& user_func) {
        local_mat.setZero();
        user_operation(e, local_mat);
        
    };


    for(const Element* e : data.elements)
    {
        //int order = e->get_geometry_order() + space_1->get_basis_order() + space_2->get_basis_order();
        Basis_Shape b = to_basis_shape(e->get_geometry());
        const FEM_Space* shape_space_1 = data.space_1->get_basis_space(b);
        const FEM_Space* shape_space_2 = data.space_2->get_basis_space(b);

        size_t col_size = shape_space_1->get_n_dof();
        size_t row_size = shape_space_1->get_n_dof();


        // dispatch to fixed-size, write logic once
        if      (col_size == 3 && row_size == 3) { Matrix<3,3> m; assemble_local(e, m); }
        else if (col_size == 4 && row_size == 4) { Matrix<4,4> m; assemble_local(e, m); }
        else if (col_size == 4 && row_size == 6) { Matrix<4,6> m; assemble_local(e, m); }
        else if (col_size == 6 && row_size == 4) { Matrix<6,4> m; assemble_local(e, m); }
        else if (col_size == 6 && row_size == 6) { Matrix<6,6> m; assemble_local(e, m); }
        else                                     { MatrixXd    m(col_size, row_size); assemble_local(e, m); }
    }



    return true;
}





}