#include "math/operator/interpolator.h"

#include "math/fem/element_data.h"


using namespace simu;

/*
template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator__H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
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
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator__H1_to_Hcurl, interpolate_element)
*/







template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator__H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;
    constexpr int S = (C == Eigen::Dynamic) ? Eigen::Dynamic : C / phy_dim;


    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTERPOLATOR_ID], [&]{ do_once(e_data.space_1, e_data.space_2, e_data.dof_manager, e_data.entity_size); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* target_space = e_data.shape_space_1;  // Hcurl space
        const FEM_Space* source_space = e_data.shape_space_2;  // H1    space

        int n_dof_per_edge = target_space->get_n_dof_per_edge();
        int n_dof_per_face = target_space->get_n_dof_per_face();
        int n_dof_per_cell = target_space->get_n_dof_per_cell();

        const size_t * e_node = e->get_node_idx();
        Basis_Shape e_shape = e_data.b_shape;

        
        auto process_entity = [&](
            const std::vector<Ref_Coord>&            coords,
            const std::vector<Integration_Point>&    i_r_rule,
            int                                      entity_dim,
            int                                      entity_idx, 
            int                                      n_dof, 
            size_t                                   row_offset, 
            auto&&                                   compute_t_phys)    // (J, dof_info_row, coord) -> Vector<phy_dim>
        {
            Matrix<Eigen::Dynamic, ref_dim> dof_info(coords.size()*n_dof, ref_dim);
            target_space->dof_signature(entity_dim, entity_idx, coords, dof_info);

            std::vector<Vector<S>>                H1_k(coords.size(), Vector<S>(e_data.cols/phy_dim));
            std::vector<Matrix<phy_dim, ref_dim>>  J_k(coords.size());
            for (int k = 0; k < coords.size(); ++k)
            {
                //std::cout<<"ref_coord: x="<<coords[k].x<<", y="<<coords[k].y<<", z="<<coords[k].z<<std::endl;
                source_space->get_basis_s(coords[k], H1_k[k]);
                J_k[k] = e_data.get_J(coords[k]);
            }
            //std::cout<<"coord size: "<<coords.size()<<std::endl;

            for (int j = 0; j < n_dof; ++j)
            {
                const size_t row = row_offset + j;

                for (int k = 0; k < coords.size(); ++k)
                {
                    Vector<phy_dim> t_phys = compute_t_phys(J_k[k], dof_info.row(j*coords.size()+k), coords[k]);
                    double w = i_r_rule[k].weight;

                    for (Eigen::Index m = 0; m < H1_k[k].size(); ++m)
                        for (int d = 0; d < phy_dim; ++d)
                            //element_matrix(row, m * phy_dim + d) += w * H1_k[k][m] * t_phys[d];
                            element_matrix(row, d*H1_k[k].size() + m) += w * H1_k[k][m] * t_phys[d];
                }
            }
        };


        // ---- edges ----
        if(n_dof_per_edge>0){
            const std::vector<Integration_Point>& i_r_edge = Integration::get_rule(get_edge(e_shape), target_space->get_basis_order()+e->get_geometry_order());
            for (int i = 0; i < e->get_n_edge(); ++i) 
            {
                // skip edge if previously visited.
                if(e_data.dof_manager->get_edge_id(e_shape, e_node, i).exist) continue;
                auto coords = e->edge_map(i_r_edge, i);
                process_entity(coords, i_r_edge, 1, i, n_dof_per_edge, i*n_dof_per_edge,
                    [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                        return J * t.transpose();
                    }
                );
                e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_edge(e_data.b_shape, e_node, i); });
                
            }
        }

        // ---- faces (for 3D elements only) ----
        if constexpr (ref_dim == 3)
        {
            if(n_dof_per_face>0){
                const size_t offset_face = n_dof_per_edge * e->get_n_edge();
                const std::vector<Integration_Point>& i_r_face = Integration::get_rule(get_face(e_shape), target_space->get_basis_order()+e->get_geometry_order());
                for (int i = 0; i < e->get_n_face(); ++i)
                {
                    // skip face if previously visited.
                    if(e_data.dof_manager->get_face_id(e_shape, e_node, i).exist) continue;
                    auto coords = e->face_map(i_r_face, i);
                    Vector<ref_dim> e1, e2;
                    e->face_ref_edge(i, e1, e2);
                    process_entity(coords, i_r_face, 2, i, n_dof_per_face, offset_face+i*n_dof_per_face,
                        [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                            Vector<phy_dim> Jt    = J * t.transpose();
                            const double    scale = ((J*e1).cross(J*e2)).norm();
                            return Jt / Jt.norm() * scale;
                        }
                    );
                    e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_face(e_data.b_shape, e_node, i); });
                }
            }
        }

        // ---- cell ----
        if(n_dof_per_cell>0){
            const size_t offset_cell = n_dof_per_edge * e->get_n_edge() + n_dof_per_face * e->get_n_face();
            const std::vector<Integration_Point>& i_r_cell = Integration::get_rule(e_shape, target_space->get_basis_order()+e->get_geometry_order());
            std::vector<Ref_Coord> coords;
            coords.reserve(i_r_cell.size());
            for (const auto& ip : i_r_cell) coords.push_back(ip.coord);

            process_entity(coords, i_r_cell, ref_dim, 0, n_dof_per_cell, offset_cell,
                [&](auto const& J, auto const& di, auto const& c) -> Vector<phy_dim> {
                    Vector<phy_dim> Jt        = J * di.transpose();
                    const double    abs_det_J = std::abs(e_data.get_det_J(c));
                    return Jt / Jt.norm() * abs_det_J;
                }
            );
        }
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator__H1_to_Hcurl, interpolate_element)










template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator__grad_H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTERPOLATOR_ID], [&]{ do_once(e_data.space_1, e_data.space_2, e_data.dof_manager, e_data.entity_size); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* target_space = e_data.shape_space_1;  // Hcurl space
        const FEM_Space* source_space = e_data.shape_space_2;  // H1    space

        int n_dof_per_edge = target_space->get_n_dof_per_edge();
        int n_dof_per_face = target_space->get_n_dof_per_face();
        int n_dof_per_cell = target_space->get_n_dof_per_cell();

        const size_t * e_node = e->get_node_idx();
        Basis_Shape e_shape = e_data.b_shape;

        
        auto process_entity = [&](
            const std::vector<Ref_Coord>&            coords,
            const std::vector<Integration_Point>&    i_r_rule,
            int                                      entity_dim,
            int                                      entity_idx, 
            int                                      n_dof, 
            size_t                                   row_offset, 
            auto&&                                   compute_t_phys)    // (J, dof_info_row, coord) -> Vector<phy_dim>
        {
            Matrix<Eigen::Dynamic, ref_dim> dof_info(coords.size()*n_dof, ref_dim);
            target_space->dof_signature(entity_dim, entity_idx, coords, dof_info);

            std::vector<Matrix<C,phy_dim>>   grad_H1_k(coords.size(), Matrix<C,ref_dim>(e_data.cols, ref_dim));
            std::vector<Matrix<phy_dim, ref_dim>>  J_k(coords.size());
            for (int k = 0; k < coords.size(); ++k)
            {
                J_k[k] = e_data.get_J(coords[k]);
                Matrix<ref_dim, phy_dim> inv_J = e_data.get_inv_J(coords[k]);
                Matrix<C,ref_dim> grad_H1;
                source_space->get_ED_basis_v(coords[k], grad_H1);
                grad_H1_k[k] = grad_H1 * inv_J;
            }

            for (int j = 0; j < n_dof; ++j)
            {
                const size_t row = row_offset + j;

                for (int k = 0; k < coords.size(); ++k)
                {
                    Vector<phy_dim> t_phys = compute_t_phys(J_k[k], dof_info.row(j*coords.size()+k), coords[k]);
                    double w = i_r_rule[k].weight;

                    element_matrix.row(row) += w * (grad_H1_k[k] * t_phys).transpose();
                }
            }
        };


        // ---- edges ----
        if(n_dof_per_edge>0){
            const std::vector<Integration_Point>& i_r_edge = Integration::get_rule(get_edge(e_shape), target_space->get_basis_order()+e->get_geometry_order());
            for (int i = 0; i < e->get_n_edge(); ++i) 
            {
                // skip edge if previously visited.
                if(e_data.dof_manager->get_edge_id(e_shape, e_node, i).exist) continue;
                auto coords = e->edge_map(i_r_edge, i);
                process_entity(coords, i_r_edge, 1, i, n_dof_per_edge, i*n_dof_per_edge,
                    [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                        return J * t.transpose();
                    }
                );
                e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_edge(e_data.b_shape, e_node, i); });
                
            }
        }

        // ---- faces (for 3D elements only) ----
        if constexpr (ref_dim == 3)
        {
            if(n_dof_per_face>0){
                const size_t offset_face = n_dof_per_edge * e->get_n_edge();
                const std::vector<Integration_Point>& i_r_face = Integration::get_rule(get_face(e_shape), target_space->get_basis_order()+e->get_geometry_order());
                for (int i = 0; i < e->get_n_face(); ++i)
                {
                    // skip face if previously visited.
                    if(e_data.dof_manager->get_face_id(e_shape, e_node, i).exist) continue;
                    auto coords = e->face_map(i_r_face, i);
                    Vector<ref_dim> e1, e2;
                    e->face_ref_edge(i, e1, e2);
                    process_entity(coords, i_r_face, 2, i, n_dof_per_face, offset_face+i*n_dof_per_face,
                        [&](auto const& J, auto const& t, auto const&) -> Vector<phy_dim> {
                            Vector<phy_dim> Jt    = J * t.transpose();
                            const double    scale = ((J*e1).cross(J*e2)).norm();
                            return Jt / Jt.norm() * scale;
                        }
                    );
                    e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_face(e_data.b_shape, e_node, i); });
                }
            }
        }

        // ---- cell ----
        if(n_dof_per_cell>0){
            const size_t offset_cell = n_dof_per_edge * e->get_n_edge() + n_dof_per_face * e->get_n_face();
            const std::vector<Integration_Point>& i_r_cell = Integration::get_rule(e_shape, target_space->get_basis_order()+e->get_geometry_order());
            std::vector<Ref_Coord> coords;
            coords.reserve(i_r_cell.size());
            for (const auto& ip : i_r_cell) coords.push_back(ip.coord);

            process_entity(coords, i_r_cell, ref_dim, 0, n_dof_per_cell, offset_cell,
                [&](auto const& J, auto const& di, auto const& c) -> Vector<phy_dim> {
                    Vector<phy_dim> Jt        = J * di.transpose();
                    const double    abs_det_J = std::abs(e_data.get_det_J(c));
                    return Jt / Jt.norm() * abs_det_J;
                }
            );
        }
        
        
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator__grad_H1_to_Hcurl, interpolate_element)











template<int phy_dim, int ref_dim, typename Mat_Type>
void Identity_Mapping::direct_mapping(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    std::call_once(e_data.check->interpolator_check[INTERPOLATOR_ID], [&]{ do_once(e_data.space_1, e_data.space_2, e_data.dof_manager, e_data.entity_size); }); 

    const Element* e = e_data.e;
    const size_t * e_node = e->get_node_idx();
    Basis_Shape e_shape = e_data.b_shape;

    const FEM_Space* fe_space = e_data.shape_space_1;

    int n_dof_per_node = fe_space->get_n_dof_per_node();
    int n_dof_per_edge = fe_space->get_n_dof_per_edge();
    int n_dof_per_face = fe_space->get_n_dof_per_face();
    int n_dof_per_cell = fe_space->get_n_dof_per_cell();

    int local_dof_idx = 0;

    if(n_dof_per_node > 0) for(int i=0; i<e->get_n_node(); ++i) 
    {
        // skip node if previously visited.
        if(e_data.dof_manager->get_node_id(e_shape, e_node, i).exist){ local_dof_idx+=n_dof_per_node; continue; }
        for(int j=0; j<n_dof_per_node; ++j){
            element_matrix(local_dof_idx, local_dof_idx) += 1;
            local_dof_idx++;
        }
        e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_node(e_data.b_shape, e_node, i); });
    }

    if(n_dof_per_edge > 0) for(int i=0; i<e->get_n_edge(); ++i) 
    {
        // skip edge if previously visited.
        if(e_data.dof_manager->get_edge_id(e_shape, e_node, i).exist) { local_dof_idx+=n_dof_per_edge; continue; }
        for(int j=0; j<n_dof_per_edge; ++j){
            element_matrix(local_dof_idx, local_dof_idx) += 1;
            local_dof_idx++;
        }
        e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_edge(e_data.b_shape, e_node, i); });
    }

    if(n_dof_per_face > 0) for(int i=0; i<e->get_n_face(); ++i) 
    {
        // skip face if previously visited.
        if(e_data.dof_manager->get_face_id(e_shape, e_node, i).exist) { local_dof_idx+=n_dof_per_face; continue; }
        for(int j=0; j<n_dof_per_face; ++j){
            element_matrix(local_dof_idx, local_dof_idx) += 1;
            local_dof_idx++;
        }
        e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_face(e_data.b_shape, e_node, i); });
    }

    for(int i=0; i<e->get_n_cell(); ++i) 
    {
        // skip cell if previously visited.
        if(e_data.dof_manager->get_cell_id(e_shape, e_node, i).exist) { local_dof_idx+=n_dof_per_cell; continue; }
        for(int j=0; j<n_dof_per_cell; ++j){
            element_matrix(local_dof_idx, local_dof_idx) += 1;
            local_dof_idx++;
        }
        e_data.dof_manager->pending_operation([&e_data, e_node, i]() { e_data.dof_manager->register_cell(e_data.b_shape, e_node, i); });
    }
    


}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Identity_Mapping, direct_mapping)
