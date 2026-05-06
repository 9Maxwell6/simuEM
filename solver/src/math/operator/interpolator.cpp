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
        const FEM_Space* source_space = e_data.shape_space_1;  // H1 space      N dof
        const FEM_Space* target_space = e_data.shape_space_2;  // Hcurl space   M dof

        

        // edge integral
        const std::vector<Integration_Point>& i_r_edge = Integration::get_rule(get_edge(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
        int n_dof_per_edge = target_space->get_n_dof_per_edge();


        Vector<R> H1_basis(e_data.rows);


        for (int i=0; i<e->get_n_edge(); ++i) 
        {
            
            std::vector<Ref_Coord> edge_coord = e->edge_map(i_r_edge, i);    // map coord on 1D unit edge to i-th edge on reference tetrahedron
            Matrix<Eigen::Dynamic, ref_dim> dof_info(edge_coord.size(),ref_dim);

            for(int j=0; j<n_dof_per_edge; ++j)
            {
                int row_idx = i*n_dof_per_edge +j;
                target_space->dof_signature(1,j,edge_coord, dof_info);  // edge tangent
 
                for(int k=0; k<edge_coord.size(); ++k)
                {
                    source_space->get_basis_s(edge_coord[k], H1_basis);
                    const Matrix<phy_dim, ref_dim>& J = e_data.get_J(edge_coord[k]);   
                    Vector<phy_dim> t_phys = J*dof_info.row(k).transpose();

                    for (Eigen::Index m = 0; m < H1_basis.size(); ++m)
                    {
                        for(int d=0; d<phy_dim; ++d)
                        {
                            double tangential = H1_basis[m] * t_phys[d];

                            int col_idx = m*phy_dim + d;

                            element_matrix(row_idx, col_idx) += i_r_edge[k].weight * tangential;
                        }
                    }
                }
            }
        }

        size_t row_offset = n_dof_per_edge*e->get_n_edge();

        // face integral
        const std::vector<Integration_Point>& i_r_face = Integration::get_rule(get_face(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
        int n_dof_per_face = target_space->get_n_dof_per_face();


        for (int i=0; i<e->get_n_face(); ++i) 
        {
            std::vector<Ref_Coord> face_coord = e->face_map(i_r_face, i);    // map coord on 2D reference triangle to i-th face on reference tetrahedron
            Matrix<Eigen::Dynamic, ref_dim> dof_info(face_coord.size(),ref_dim);

            for(int j=0; j<n_dof_per_edge; ++j)
            {
                int row_idx = row_offset + i*n_dof_per_edge +j;
                target_space->dof_signature(2,j,face_coord, dof_info);  // face tangent
 
                for(int k=0; k<face_coord.size(); ++k)
                {
                    source_space->get_basis_s(face_coord[k], H1_basis);
                    const Matrix<phy_dim, ref_dim>& J = e_data.get_J(face_coord[k]);   
                    Vector<phy_dim> t_phys = J*dof_info.row(k).transpose();

                    for (Eigen::Index m = 0; m < H1_basis.size(); ++m)
                    {
                        for(int d=0; d<phy_dim; ++d)
                        {
                            double tangential = H1_basis[m] * t_phys[d];

                            int col_idx = m*phy_dim + d;

                            element_matrix(row_idx, col_idx) += i_r_face[k].weight * tangential;
                        }
                    }
                }
            }
        }


        // cell integral
        const std::vector<Integration_Point>& i_r_cell = Integration::get_rule(e_data.b_shape, target_space->get_basis_order()+e->get_geometry_order());
        int n_dof_per_cell = target_space->get_n_dof_per_cell();

        if constexpr(ref_dim == phy_dim){

            for (int i=0; i<e->get_n_face(); ++i) 
            {
                std::vector<Ref_Coord> face_coord = e->face_map(i_r_cell, i);    // map coord on 2D reference triangle to i-th face on reference tetrahedron
                Matrix<Eigen::Dynamic, ref_dim> dof_info(face_coord.size(),ref_dim);

                for(int j=0; j<n_dof_per_edge; ++j)
                {
                    int row_idx = row_offset + i*n_dof_per_edge +j;
                    target_space->dof_signature(3,j,face_coord, dof_info);  
    
                    for(int k=0; k<face_coord.size(); ++k)
                    {
                        source_space->get_basis_s(face_coord[k], H1_basis);
                        double abs_det_J = std::abs(e_data.get_det_J(face_coord[k]));   
                        Vector<phy_dim> t_phys = abs_det_J*dof_info.row(k).transpose();

                        for (Eigen::Index m = 0; m < H1_basis.size(); ++m)
                        {
                            for(int d=0; d<phy_dim; ++d)
                            {
                                double tangential = H1_basis[m] * t_phys[d];

                                int col_idx = m*phy_dim + d;

                                element_matrix(row_idx, col_idx) += i_r_cell[k].weight * tangential;
                            }
                        }
                    }
                }
            }
        }



    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator_H1_to_Hcurl, interpolate_element)


