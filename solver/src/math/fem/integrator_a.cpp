#include "math/fem/integrator_a.h"

#include "math/fem/element_data.h"



using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_S__S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Vector<R> basis(e_data.rows);

        for(const Integration_Point& i_p : i_p_list)
        {
            trial_space->get_basis_s(i_p.coord, basis);
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));

            element_matrix += coeff * i_p.weight * abs_det_J * basis * basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_S__S, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_grad_S__grad_S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
    
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order()-1 + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
    
        Matrix<R,ref_dim> grad_basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_grad_basis(e_data.rows, phy_dim);

        for(const Integration_Point& i_p : i_p_list)
        {            
            trial_space->get_ED_basis_v(i_p.coord, grad_basis);
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& inv_J = e_data.get_inv_J(i_p.coord);

            phy_grad_basis = grad_basis*inv_J;

            element_matrix += coeff *i_p.weight * abs_det_J * phy_grad_basis * phy_grad_basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_grad_S__grad_S, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator_H1__s_V__grad_S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 

        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);

        if constexpr ((R == Eigen::Dynamic && C == Eigen::Dynamic) || (phy_dim*R == C)) {
            
            if(phy_dim*e_data.rows != e_data.cols) { 
                Logger::error("Integrator_H1__s_V__grad_S::assemble_element_matrix: incorrect matrix dimensions.");
                return;
            }
            
            Vector<R> basis;

            Matrix<R, ref_dim> grad_basis(e_data.rows, ref_dim);
            Matrix<R, phy_dim> phy_grad_basis(e_data.rows, phy_dim);

            for(const Integration_Point& i_p : i_p_list)
            {            
                test__space->get_basis_s(i_p.coord, basis);
                trial_space->get_ED_basis_v(i_p.coord, grad_basis);

                
                double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
                const MatrixXd& inv_J = e_data.get_inv_J(i_p.coord);

                phy_grad_basis = grad_basis*inv_J;

                for(size_t i=0; i<phy_dim; ++i)
                {
                    element_matrix.block(0, i*e_data.rows, e_data.rows, e_data.rows) += coeff * i_p.weight * abs_det_J * basis * phy_grad_basis.col(i).transpose();
                }
            }
        }
    
        // TODO:  need to make sure the dof index from the r.h.s linear form match the dof from l.h.s bilinear form.

        if(coeff!=0){
            std::cout<<element_matrix<<std::endl;
            std::cout<<"====================================="<<std::endl;
            const size_t * e_list = e->get_node_idx();
            const Node& n0 = e_data.mesh->get_node(e_list[0]);
            const Node& n1 = e_data.mesh->get_node(e_list[1]);
            const Node& n2 = e_data.mesh->get_node(e_list[2]);
            const Node& n3 = e_data.mesh->get_node(e_list[3]);
            std::cout<<"x="<<n0.x<<", y="<<n0.y<<", z="<<n0.z<<std::endl;
            std::cout<<"x="<<n1.x<<", y="<<n1.y<<", z="<<n1.z<<std::endl;
            std::cout<<"x="<<n2.x<<", y="<<n2.y<<", z="<<n2.z<<std::endl;
            std::cout<<"x="<<n3.x<<", y="<<n3.y<<", z="<<n3.z<<std::endl;
            exit(0);
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator_H1__s_V__grad_S, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_curl_V__curl_V::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (((phy_dim == 3 && (ref_dim == 3 || ref_dim ==2)) || (phy_dim == 2 && (phy_dim == ref_dim))) && R == C) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 

        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order()-1 + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Matrix<R, ref_dim> curl_basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_curl_basis(e_data.rows, phy_dim);

        for(const Integration_Point& i_p : i_p_list)
        {            
            trial_space->get_ED_basis_v(i_p.coord, curl_basis);
            
            double det_J = e_data.get_det_J(i_p.coord);
            const Matrix<phy_dim, ref_dim>& J = e_data.get_J(i_p.coord);            

            if constexpr  (phy_dim == 3){
                phy_curl_basis = curl_basis * J.transpose()/det_J;
            }else if constexpr (phy_dim == ref_dim && ref_dim == 2){
                phy_curl_basis = curl_basis/det_J;   // TODO: need verify
            }else{
                Logger::error("Integrator__s_curl_V__curl_V::assemble_element_matrix: reference element dimension = 1 not available.");
                return;
            }
            element_matrix += coeff * i_p.weight * std::abs(det_J) * phy_curl_basis * phy_curl_basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_curl_V__curl_V, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_V__V::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 

        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order();
        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Matrix<R, ref_dim> basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_basis(e_data.rows, phy_dim);
        for(const Integration_Point& i_p : i_p_list)
        {
            trial_space->get_basis_v(i_p.coord, basis);
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& J_inv = e_data.get_inv_J(i_p.coord);   

            phy_basis = basis * J_inv;

            element_matrix += coeff * i_p.weight * abs_det_J * phy_basis * phy_basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_V__V, assemble_element_matrix, double)






template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_V__grad_S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;
    
    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 

        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order()-1 + test__space->get_basis_order();
        
        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Matrix<R, ref_dim> domain_grad_basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_domain_grad_basis(e_data.rows, phy_dim);


        Matrix<C, ref_dim> range_basis(e_data.cols, ref_dim);
        Matrix<C, phy_dim> phy_range_basis(e_data.cols, phy_dim);
        for(const Integration_Point& i_p : i_p_list)
        {
            test__space->get_ED_basis_v(i_p.coord, domain_grad_basis);
            trial_space->get_basis_v(i_p.coord, range_basis);
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& J_inv = e_data.get_inv_J(i_p.coord);   

            phy_domain_grad_basis = domain_grad_basis * J_inv;
            phy_range_basis = range_basis * J_inv;

            element_matrix += coeff * i_p.weight * abs_det_J * phy_domain_grad_basis * phy_range_basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_V__grad_S, assemble_element_matrix, double)





template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_grad_S__V::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 

        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);

        Matrix<R, ref_dim> domain_basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_domain_basis(e_data.rows, phy_dim);

        

        Matrix<C, ref_dim> range_grad_basis(e_data.cols, ref_dim);
        Matrix<C, phy_dim> phy_range_grad_basis(e_data.cols, phy_dim);

        for(const Integration_Point& i_p : i_p_list)
        {
            test__space->get_basis_v(i_p.coord, domain_basis);         // row
            trial_space->get_ED_basis_v(i_p.coord, range_grad_basis);  // column
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& J_inv = e_data.get_inv_J(i_p.coord);   

            phy_domain_basis = domain_basis * J_inv;
            phy_range_grad_basis = range_grad_basis * J_inv;

            element_matrix += coeff * i_p.weight * abs_det_J * phy_domain_basis * phy_range_grad_basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Integrator__s_grad_S__V, assemble_element_matrix, double)



