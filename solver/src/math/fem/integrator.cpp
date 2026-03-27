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

        if(e_data.i_r_list->size() < order) e_data.i_r_list->resize(order+1, nullptr);
        
        const std::vector<Integration_Point>* i_p_list = (*e_data.i_r_list)[order];
        if(!i_p_list) {
            i_p_list = &Integration::get_integration_points(e_data.b_shape, order);
            (*e_data.i_r_list)[order] = i_p_list;
        }

        
        
        for(int i=0; i<i_p_list->size(); ++i)
        {
            Vector<R> basis;
            trial_space->get_basis_s((*i_p_list)[i], basis);
            
            double det_J = e_data.get_det_Jacobian((*i_p_list)[i]);

            if(det_J==0){
                const size_t * e_list = e->get_node_idx();
                std::cout<<e_list[0]<<std::endl;
                std::cout<<e_list[1]<<std::endl;
                std::cout<<e_list[2]<<std::endl;
                std::cout<<e_list[3]<<std::endl;
                std::cout<<"-----------------------------------"<<std::endl;

                const Node& n0 = e_data.mesh->get_node(e_list[0]);
                const Node& n1 = e_data.mesh->get_node(e_list[1]);
                const Node& n2 = e_data.mesh->get_node(e_list[2]);
                const Node& n3 = e_data.mesh->get_node(e_list[3]);

                std::cout<<"x: "<<n0.x<<",  y: "<<n0.y<<",  z: "<<n0.z<<std::endl;
                std::cout<<"x: "<<n1.x<<",  y: "<<n1.y<<",  z: "<<n1.z<<std::endl;
                std::cout<<"x: "<<n2.x<<",  y: "<<n2.y<<",  z: "<<n2.z<<std::endl;
                std::cout<<"x: "<<n3.x<<",  y: "<<n3.y<<",  z: "<<n3.z<<std::endl;
            }
            
        }

        

        
    }
}

INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_S__S, assemble_element_matrix, double)


