#include "math/operator/interpolator.h"

#include "math/fem/element_data.h"


using namespace simu;

/*
template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator_H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTERPOLATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* source_space = e_data.shape_space_1;  // H1 space      N dof
        const FEM_Space* target_space = e_data.shape_space_2;  // Hcurl space   M dof

        int n_dof_per_edge = target_space->get_n_dof_per_edge();
        int n_dof_per_face = target_space->get_n_dof_per_face();
        int n_dof_per_cell = target_space->get_n_dof_per_cell();

        // edge integral
        const std::vector<Integration_Point>& i_r_edge = Integration::get_rule(get_edge(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
        


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
                            element_matrix(row_idx, m*phy_dim + d) += i_r_edge[k].weight * H1_basis[m] * t_phys[d];
                        }
                    }
                }
            }
        }

        size_t row_offset = n_dof_per_edge*e->get_n_edge();

        // only 3D elements need face integral
        if constexpr(ref_dim == 3){

            // face integral
            const std::vector<Integration_Point>& i_r_face = Integration::get_rule(get_face(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
            

            Vector<ref_dim> ref_edge_1;
            Vector<ref_dim> ref_edge_2;

            for (int i=0; i<e->get_n_face(); ++i) 
            {
                std::vector<Ref_Coord> face_coord = e->face_map(i_r_face, i);    // map coord on 2D reference triangle to i-th face on reference tetrahedron
                Matrix<Eigen::Dynamic, ref_dim> dof_info(face_coord.size(),ref_dim);
                e->face_ref_edge(i, ref_edge_1, ref_edge_2);

                for(int j=0; j<n_dof_per_face; ++j)
                {
                    int row_idx = row_offset + i*n_dof_per_face +j;
                    target_space->dof_signature(2,j,face_coord, dof_info);  //  face tangent
    
                    for(int k=0; k<face_coord.size(); ++k)
                    {
                        source_space->get_basis_s(face_coord[k], H1_basis);
                        const Matrix<phy_dim, ref_dim>& J = e_data.get_J(face_coord[k]);   
                        double scale                      = ((J*ref_edge_1).cross(J*ref_edge_2)).norm();
                        Vector<phy_dim> Jt                = J * dof_info.row(k).transpose();
                        Vector<phy_dim> t_phys            = Jt/Jt.norm() * scale;

                        for (Eigen::Index m = 0; m < H1_basis.size(); ++m)
                        {
                            for(int d=0; d<phy_dim; ++d)
                            {
                                element_matrix(row_idx, m*phy_dim + d) += i_r_face[k].weight * H1_basis[m] * t_phys[d];
                            }
                        }
                    }
                }
            }
        }

        row_offset = n_dof_per_edge*e->get_n_edge() + n_dof_per_face*e->get_n_face();


        // cell integral
        const std::vector<Integration_Point>& i_r_cell = Integration::get_rule(e_data.b_shape, target_space->get_basis_order()+e->get_geometry_order());
        std::vector<Ref_Coord> cell_coord;
        cell_coord.reserve(i_r_cell.size());
        for (const auto& ip : i_r_cell) cell_coord.push_back(ip.coord);

        Matrix<Eigen::Dynamic, ref_dim> dof_info(cell_coord.size(),ref_dim);

        // each element only has 1 cell
        for(int j=0; j<n_dof_per_cell; ++j)
        {
            int row_idx = row_offset + j;
            target_space->dof_signature(ref_dim ,j,cell_coord, dof_info);  

            for(int k=0; k<cell_coord.size(); ++k)
            {
                source_space->get_basis_s(cell_coord[k], H1_basis);
                const Matrix<phy_dim, ref_dim>& J = e_data.get_J(cell_coord[k]);   
                double abs_det_J                  = std::abs(e_data.get_det_J(cell_coord[k]));   
                Vector<phy_dim> Jt                = J * dof_info.row(k).transpose();
                Vector<phy_dim> t_phys            = (Jt / Jt.norm()) * abs_det_J;

                for (Eigen::Index m = 0; m < H1_basis.size(); ++m)
                {
                    for(int d=0; d<phy_dim; ++d)
                    {
                        element_matrix(row_idx, m*phy_dim + d) += i_r_cell[k].weight * H1_basis[m] * t_phys[d];
                    }
                }
            }
        }
        
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator_H1_to_Hcurl, interpolate_element)
*/







template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator_H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTERPOLATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* source_space = e_data.shape_space_1;  // H1 space      N dof
        const FEM_Space* target_space = e_data.shape_space_2;  // Hcurl space   M dof

        int n_dof_per_edge = target_space->get_n_dof_per_edge();
        int n_dof_per_face = target_space->get_n_dof_per_face();
        int n_dof_per_cell = target_space->get_n_dof_per_cell();

        
        auto process_entity = [&](
            const std::vector<Ref_Coord>&            coords,
            const std::vector<Integration_Point>&    i_r_rule,
            int                                      dim_arg,
            int                                      n_dof, 
            size_t                                   row_offset, 
            auto&&                                   compute_t_phys)    // (J, dof_info_row, coord) -> Vector<phy_dim>
        {
            Matrix<Eigen::Dynamic, ref_dim> dof_info(coords.size(), ref_dim);

            std::vector<Vector<R>>                H1_k(coords.size(), Vector<R>(e_data.rows));
            std::vector<Matrix<phy_dim, ref_dim>>  J_k(coords.size());
            for (int k = 0; k < coords.size(); ++k)
            {
                source_space->get_basis_s(coords[k], H1_k[k]);
                J_k[k] = e_data.get_J(coords[k]);
            }

            for (int j = 0; j < n_dof; ++j)
            {
                target_space->dof_signature(dim_arg, j, coords, dof_info);
                const size_t row = row_offset + j;

                for (int k = 0; k < (int)coords.size(); ++k)
                {
                    Vector<phy_dim> t_phys = compute_t_phys(J_k[k], dof_info.row(k), coords[k]);
                    double w = i_r_rule[k].weight;

                    for (Eigen::Index m = 0; m < H1_k[k].size(); ++m)
                        for (int d = 0; d < phy_dim; ++d)
                            element_matrix(row, m * phy_dim + d) += w * H1_k[k][m] * t_phys[d];
                }
            }
        };


        // ---- edges ----
        const std::vector<Integration_Point>& i_r_edge = Integration::get_rule(get_edge(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
        for (int i = 0; i < e->get_n_edge(); ++i) 
        {
            auto coords = e->edge_map(i_r_edge, i);
            process_entity(coords, i_r_edge, 1, n_dof_per_edge, i*n_dof_per_edge,
                [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                    return J * t.transpose();
                }
            );
        }

        // ---- faces (for 3D elements only) ----
        const size_t off_face = n_dof_per_edge * e->get_n_edge();
        if constexpr (ref_dim == 3)
        {
            const std::vector<Integration_Point>& i_r_face = Integration::get_rule(get_face(e_data.b_shape), target_space->get_basis_order()+e->get_geometry_order());
            for (int i = 0; i < e->get_n_face(); ++i)
            {
                auto coords = e->face_map(i_r_face, i);
                Vector<ref_dim> e1, e2;
                e->face_ref_edge(i, e1, e2);
                process_entity(coords, i_r_face, 2, n_dof_per_face, off_face+i*n_dof_per_face,
                    [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                        Vector<phy_dim> Jt    = J * t.transpose();
                        const double    scale = ((J*e1).cross(J*e2)).norm();
                        return Jt / Jt.norm() * scale;
                    }
                );
            }
        }

        // ---- cell ----
        const size_t off_cell = off_face + n_dof_per_face * e->get_n_face();
        const std::vector<Integration_Point>& i_r_cell = Integration::get_rule(e_data.b_shape, target_space->get_basis_order()+e->get_geometry_order());
        std::vector<Ref_Coord> coords;
        coords.reserve(i_r_cell.size());
        for (const auto& ip : i_r_cell) coords.push_back(ip.coord);

        process_entity(coords, i_r_cell, ref_dim, n_dof_per_cell, off_cell,
            [&](auto const& J, auto const& di, auto const& c) -> Vector<phy_dim> {
                Vector<phy_dim> Jt        = J * di.transpose();
                const double    abs_det_J = std::abs(e_data.get_det_J(c));
                return Jt / Jt.norm() * abs_det_J;
            }
        );
        
        
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator_H1_to_Hcurl, interpolate_element)


