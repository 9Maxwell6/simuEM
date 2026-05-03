#include "math/operator/interpolator.h"

#include "math/fem/element_data.h"


using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator_H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* source_space = e_data.shape_space_1;  // H1 space
        const FEM_Space* target_space = e_data.shape_space_2;  // Hcurl space

        // NxM dof
        int n_dof_1 = source_space->get_n_dof();  // N dof
        int n_dof_2 = target_space->get_n_dof();  // M dof


        // edge integral
        const std::vector<Integration_Point>& i_r = Integration::get_rule(get_edge(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
        int n_dof_per_edge = target_space->get_n_dof_per_edge();


        Vector<R> H1_basis(e_data.rows);

        Matrix<Eigen::Dynamic, ref_dim> all_edge_tangent(e->get_n_edge(), e->get_dim());
        e->tangent(all_edge_tangent); // get tangent vector for all edge of reference element. the magnitude is length of edge.


        for(const Integration_Point& i_p : i_r)
        {
            std::vector<Ref_Coord> all_edge_coord =  e->edge_map(i_p.coord);

            for (int i=0; i<all_edge_coord.size(); ++i) 
            {
                for(int j=0; j<n_dof_per_edge; ++j)
                {
                    //TODO: implement edge_poly per FEM_space
                    //double edge_poly = target_space->edge_poly(i_p.coord, j);
                    int row_idx = i*n_dof_per_edge +j;
                    Ref_Coord& edge_coord = all_edge_coord[i];
                    auto edge_tangent = all_edge_tangent.row(i).transpose();

                    source_space->get_basis_s(edge_coord, H1_basis);
                    const Matrix<phy_dim, ref_dim>& J = e_data.get_J(edge_coord);   
                    Vector<phy_dim> t_phys = J*edge_tangent; 

                    for (Eigen::Index k = 0; k < H1_basis.size(); ++k)
                    {
                        for(int d=0; d<phy_dim; ++d)
                        {
                            double tangential = H1_basis[k] * t_phys[d];

                            int col_idx = k*phy_dim + d;

                            element_matrix(row_idx, col_idx) += i_p.weight * tangential;
                            //element_matrix(row_idx, col_idx) += i_p.weight * tangential * edge_poly;

                        }
                    }
                }
            }
        }



        
        
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator_H1_to_Hcurl, interpolate_element)


