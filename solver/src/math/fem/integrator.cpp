#include "math/fem/integrator.h"



using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_S__S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        const Element* e = e_data.e;
        const FEM_Space* trial_space = e_data.shape_space_1;
        const FEM_Space* test__space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        

        Vector<R> basis;
        for(const Integration_Point& i_p : i_p_list)
        {
            trial_space->get_basis_s(i_p, basis);
            
            double abs_det_J = e_data.get_abs_det_J(i_p);

            element_matrix += i_p.weight * abs_det_J * basis * basis.transpose();
        }
        element_matrix *= coeff;
    }
}
INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_S__S, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_grad_S__grad_S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        const Element* e = e_data.e;
        const FEM_Space* trial_space = e_data.shape_space_1;
        const FEM_Space* test__space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order()-1 + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
    
        Matrix<R,ref_dim> grad_basis;
        Matrix<R, phy_dim> phy_grad_basis;
        for(const Integration_Point& i_p : i_p_list)
        {            
            trial_space->get_ED_basis_v(i_p, grad_basis);
            
            double abs_det_J = e_data.get_abs_det_J(i_p);
            const Matrix<ref_dim, phy_dim>& inv_J = e_data.get_inv_J(i_p);

            phy_grad_basis = grad_basis*inv_J;

            element_matrix += i_p.weight * abs_det_J * phy_grad_basis * phy_grad_basis.transpose();
        }
        element_matrix *= coeff;
    }
}
INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_grad_S__grad_S, assemble_element_matrix, double)




template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_V__grad_S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        const Element* e = e_data.e;
        const FEM_Space* trial_space = e_data.shape_space_1;
        const FEM_Space* test__space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order()-1;

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);

        if constexpr (R == Eigen::Dynamic && C ==Eigen::Dynamic) {
            if(phy_dim*element_matrix.rows()!=element_matrix.cols()) {return;}  // error
            VectorXd basis;

            MatrixXd grad_basis;
            MatrixXd phy_grad_basis;
            
            for(const Integration_Point& i_p : i_p_list)
            {            
                trial_space->get_basis_s(i_p, basis);
                test__space->get_ED_basis_v(i_p, grad_basis);

                
                double abs_det_J = e_data.get_abs_det_J(i_p);
                const MatrixXd& inv_J = e_data.get_inv_J(i_p);

                phy_grad_basis = grad_basis*inv_J;

                element_matrix += i_p.weight * abs_det_J * phy_grad_basis * phy_grad_basis.transpose();
            }
            element_matrix *= coeff;
        }else if constexpr (phy_dim*R == C) {
            Vector<R> basis;

            Matrix<R, ref_dim> grad_basis;
            Matrix<R, phy_dim> phy_grad_basis;
            
            for(const Integration_Point& i_p : i_p_list)
            {            
                trial_space->get_basis_s(i_p, basis);
                test__space->get_ED_basis_v(i_p, grad_basis);

                
                double abs_det_J = e_data.get_abs_det_J(i_p);
                const Matrix<ref_dim, phy_dim>& inv_J = e_data.get_inv_J(i_p);

                phy_grad_basis = grad_basis*inv_J;

                element_matrix += i_p.weight * abs_det_J * phy_grad_basis * phy_grad_basis.transpose();
            }
            element_matrix *= coeff;

        }

        

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
INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_V__grad_S, assemble_element_matrix, double)



