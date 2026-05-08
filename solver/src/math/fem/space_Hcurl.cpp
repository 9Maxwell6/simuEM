#include "math/fem/space_Hcurl.h"

#include <stdexcept>


using namespace simu;


Hcurl_Space::Hcurl_Space(int dim, int p) : FEM_Space(dim, p)
{
    vdim_   = 1;
    layout_ = 0;
}


bool Hcurl_Space::add_basis_shape(Basis_Shape g)
{
    std::unique_ptr<Hcurl_Space> shape_;
    switch (g) {
        case Basis_Shape::TRIANGLE:    shape_ = std::unique_ptr<Hcurl_Space>(new Hcurl_triangle(p_));    break;
        case Basis_Shape::TETRAHEDRON: shape_ = std::unique_ptr<Hcurl_Space>(new Hcurl_tetrahedron(p_)); break;
        default: 
        {
            Logger::warning("Hcurl_Space::add_basis_shape - shape not available: return false");
            return false;
        }
    }

    // uniqueness check
    for (auto& existing : shape_Hcurl_)
    {
        if (typeid(*existing.second) == typeid(*shape_)) return false;
    }

    //shape_Hcurl_.push_back(std::move(shape_));
    shape_Hcurl_[g] = std::move(shape_);
    basis_shapes_.push_back(g);
    return true;

}



FEM_Space * Hcurl_Space::get_basis_space(Basis_Shape s) const
{
    auto it = shape_Hcurl_.find(s);
    if (it != shape_Hcurl_.end()) return it->second.get();

    Logger::error("Hcurl_tetrahedron::get_basis_space - failed: key not found, return nullptr.");
    return nullptr;
}



const std::vector<Basis_Shape>& Hcurl_Space::get_basis_shapes() const
{
    return basis_shapes_;
}


//void Hcurl_Space::get_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const { shape_Hcurl_.at(s)->get_basis_s(p, basis); }
//void Hcurl_Space::get_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const { shape_Hcurl_.at(s)->get_basis_v(p, basis); }
//void Hcurl_Space::get_ED_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const { shape_Hcurl_.at(s)->get_ED_basis_s(p, basis); }
//void Hcurl_Space::get_ED_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const { shape_Hcurl_.at(s)->get_ED_basis_v(p, basis); }




/**
 * ================================= Triangle =================================
 * 
 */

Hcurl_triangle::Hcurl_triangle(int p) : Hcurl_Space(2, p)
{
    n_node_   = 3;
    n_edge_   = 3;
    n_face_   = 0;
    n_cell_ = 1;


    n_dof_            = p*(p+2);
    n_dof_per_node_   = 0;
    n_dof_per_edge_   = p;
    n_dof_per_face_   = 0;
    n_dof_per_cell_ = p*(p-1);
}

void Hcurl_triangle::get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const 
{
    double x = coord.x;
    double y = coord.y;

    switch(p_)
    {
        case 1:
            //  with Barycentric coordinates λ
            //          λ0 = 1.0 - x - y ;
            //          λ1 = x;
            //          λ2 = y;
            // => local Whitney-1 form:  W_ij = λi*∇λj - λj*∇λi,  
            //    where i≠j are local index of triangle vertices.
            //
            if (basis.rows() != 3 || basis.cols() != 2) {
                throw std::invalid_argument("matrix must be 3x2 for p-1 H(curl) Triangle.");
            }

            basis <<  1.0-y  ,    x    ,   // Edge 0: 0 -> 1
                        y    ,  1.0-x  ,   // Edge 1: 0 -> 2
                       -y    ,    x    ;   // Edge 2: 1 -> 2
                     
            break;
        default:
            throw std::invalid_argument("Edge element not available for order:  "+std::to_string(p_));
    }

}

void Hcurl_triangle::get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    switch(p_)
    {
        case 1:
            // (2D curl operator)
            // curl [1-y  , x    ]
            // curl [y    , 1-x  ]
            // curl [-y   , x    ]
            if (basis.rows() != 3 || basis.cols() != 1) {
                throw std::invalid_argument("matrix must be 3x1 for p-1 H(curl) Triangle.");
            }

            basis << 2.0,  // Edge 0
                    -2.0,  // Edge 1
                     2.0;  // Edge 2
            break;
        default:
            throw std::invalid_argument("Edge element not available for order:  "+std::to_string(p_));
    }
    
}


void Hcurl_triangle::dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const
{
    switch (p_)
    {
    case 1:
    {
        // P is 3x3 diagonal matrix of +1/-1
        // Edge 0: 0 -> 1
        // Edge 1: 0 -> 2
        // Edge 2: 1 -> 2
        size_t idx_0 = node_idx[0];
        size_t idx_1 = node_idx[1];
        size_t idx_2 = node_idx[2];

        double s0 = (idx_0 < idx_1) ? 1. : -1.;
        double s1 = (idx_0 < idx_2) ? 1. : -1.;
        double s2 = (idx_1 < idx_2) ? 1. : -1.;

        P << s0,  0,  0,
              0, s1,  0,
              0,  0, s2;
                
        break;
    }
    default:
        Logger::warning("Triangle::compute_dof_transformation_H_curl: higher order case not available.");
        break;
    }
}



void Hcurl_triangle::dof_signature(
    int entity_dim, 
    int entity_idx, 
    const std::vector<Ref_Coord>& coord_list, 
    Eigen::Ref<MatrixXd> kernel
) const
{
    //TODO: normalization factor?
    //const double inv_sqrt2 = 1. / std::sqrt(2.);
    size_t row_idx = 0;
    switch (p_)
    {
    case 1:
        if(entity_dim == 1){   // edge
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  1.       , 0.       ;  break;   // edge 0 dof 0
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 1.       ;  break;   // edge 1 dof 1
            case 2: for(auto& s : coord_list) kernel.row(row_idx++) << -1, 1;  break;   // edge 2 dof 2
            }
        }
        break;
    case 2:
    {
        if(entity_dim == 1){   // edge
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  (1.-s.x)          ,  0.               ;   // edge 0 dof 0
                    for(auto& s : coord_list) kernel.row(row_idx++) <<  (   s.x)          ,  0.               ;   // edge 0 dof 1
                break;    
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<   0.               , (1.-s.x)          ;   // edge 1 dof 2
                    for(auto& s : coord_list) kernel.row(row_idx++) <<   0.               , (   s.x)          ;   // edge 1 dof 3
                break;
            case 2: for(auto& s : coord_list) kernel.row(row_idx++) << -(1.-s.x), (1.-s.x);   // edge 2 dof 4
                    for(auto& s : coord_list) kernel.row(row_idx++) << -(   s.x), (   s.x);   // edge 2 dof 5
                break;
            }
        }else if(entity_dim == 2){  // face
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  1., 0.; break;    // face 0 dof 6
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<  0., 1.; break;    // face 1 dof 7
            }
        }
        break;
    }
    default:
        Logger::warning("Hcurl_triangle::get_dof_signature: higher order case not available.");
    }
}


/**
 * =============================== Tetrahedron ===============================
 * 
 */

Hcurl_tetrahedron::Hcurl_tetrahedron(int p) : Hcurl_Space(3, p)
{
    n_node_   = 4;
    n_edge_   = 6;
    n_face_   = 4;
    n_cell_ = 1;

    n_dof_            = p*(p+2)*(p+3)/2;
    n_dof_per_node_   = 0;
    n_dof_per_edge_   = p;
    n_dof_per_face_   = p*(p-1);
    n_dof_per_cell_ = p*(p-1)*(p-2)/2;

}

void Hcurl_tetrahedron::get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const 
{
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;

    switch(p_)
    {
        case 1:
            // https://defelement.org/elements/examples/tetrahedron-nedelec1-lagrange-0.html
            //  with Barycentric coordinates λ
            //          λ0 = 1.0 - x - y - z;
            //          λ1 = x;
            //          λ2 = y;
            //          λ3 = z;
            // => local Whitney-1 form:  W_ij = λi*∇λj - λj*∇λi,  
            //    where i≠j are local index of tetrahedron vertices.
            // => 6x3 dense matrix
            //
            if (basis.rows() != 6 || basis.cols() != 3) {
                throw std::invalid_argument("matrix must be 6x3 for p-1 H(curl) Tetrahedron.");
            }

            basis << 1.0-y-z ,    x    ,    x    ,  // Edge 0: 0 -> 1
                        y    , 1.0-x-z ,    y    ,  // Edge 1: 0 -> 2
                        z    ,    z    , 1.0-x-y ,  // Edge 2: 0 -> 3
                       -y    ,    x    ,   0.0   ,  // Edge 3: 1 -> 2
                       -z    ,   0.0   ,    x    ,  // Edge 4: 1 -> 3
                       0.0   ,   -z    ,    y    ;  // Edge 5: 2 -> 3
            break;
        default:
            throw std::invalid_argument("Edge element not available for order:  "+std::to_string(p_));
    }

}

void Hcurl_tetrahedron::get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;
    switch(p_)
    {
        case 1:
            // curl [1-y-z, x    , x    ]
            // curl [y    , 1-x-z, y    ]
            // curl [z    , z    , 1-x-y]
            // curl [-y   , x    , 0    ]
            // curl [-z   , 0    , x    ]
            // curl [0    , -z   , y    ]
            if (basis.rows() != 6 || basis.cols() != 3) {
                throw std::invalid_argument("matrix must be 6x3 for p-1 H(curl) Tetrahedron.");
            }

            basis << 0.0, -2.0,  2.0,  // Edge 0
                     2.0,  0.0, -2.0,  // Edge 1
                    -2.0,  2.0,  0.0,  // Edge 2
                     0.0,  0.0,  2.0,  // Edge 3
                     0.0, -2.0,  0.0,  // Edge 4
                     2.0,  0.0,  0.0;  // Edge 5
            break;
        default:
            throw std::invalid_argument("Edge element not available for order:  "+std::to_string(p_));
    }
    
}


void Hcurl_tetrahedron::dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const
{
    switch (p_)
    {
    case 1:
    {
        // P is 6x6 diagonal matrix of +1/-1
        // Edge 0: 0 -> 1
        // Edge 1: 0 -> 2
        // Edge 2: 0 -> 3
        // Edge 3: 1 -> 2
        // Edge 4: 1 -> 3
        // Edge 5: 2 -> 3
        size_t idx_0 = node_idx[0];
        size_t idx_1 = node_idx[1];
        size_t idx_2 = node_idx[2];
        size_t idx_3 = node_idx[3];

        double s0 = (idx_0 < idx_1) ? 1. : -1.;
        double s1 = (idx_0 < idx_2) ? 1. : -1.;
        double s2 = (idx_0 < idx_3) ? 1. : -1.;
        double s3 = (idx_1 < idx_2) ? 1. : -1.;
        double s4 = (idx_1 < idx_3) ? 1. : -1.;
        double s5 = (idx_2 < idx_3) ? 1. : -1.;

        P << s0,  0,  0,  0,  0,  0,
              0, s1,  0,  0,  0,  0,
              0,  0, s2,  0,  0,  0,
              0,  0,  0, s3,  0,  0,
              0,  0,  0,  0, s4,  0,
              0,  0,  0,  0,  0, s5;
        break;
    }
    default:
        Logger::warning("Tetrahedron::compute_dof_transformation_H_curl: higher order case not available.");
        break;
    }
}


void Hcurl_tetrahedron::dof_signature(
    int entity_dim, 
    int entity_idx, 
    const std::vector<Ref_Coord>& coord_list, 
    Eigen::Ref<MatrixXd> kernel
) const
{
    //TODO: normalization factor?
    const double inv_sqrt2 = 1. / std::sqrt(2.);
    const double inv_sqrt3 = 1. / std::sqrt(3.);
    size_t row_idx = 0;
    switch (p_)
    {
    case 1:
        if(entity_dim == 1){   // edge
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  1.       , 0.       , 0.       ;  break;   // edge 0 dof 0
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 1.       , 0.       ;  break;   // edge 1 dof 1
            case 2: for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 0.       , 1.       ;  break;   // edge 2 dof 2
            case 3: for(auto& s : coord_list) kernel.row(row_idx++) << -1, 1, 0.       ;  break;   // edge 3 dof 3
            case 4: for(auto& s : coord_list) kernel.row(row_idx++) << -1, 0.       , 1;  break;   // edge 4 dof 4
            case 5: for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       ,-1, 1;  break;   // edge 5 dof 5
            }
        }
        break;
    case 2:
    {
        if(entity_dim == 1){   // edge
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  (1.-s.x),  0.     ,  0.     ;   // edge 0 dof 0
                    for(auto& s : coord_list) kernel.row(row_idx++) <<  (   s.x),  0.     ,  0.     ;   // edge 0 dof 1
                break;    
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<   0.     , (1.-s.x),  0.     ;   // edge 1 dof 2
                    for(auto& s : coord_list) kernel.row(row_idx++) <<   0.     , (   s.x),  0.     ;   // edge 1 dof 3
                break;
            case 2: for(auto& s : coord_list) kernel.row(row_idx++) <<   0.     ,  0.     , (1.-s.x);   // edge 2 dof 4
                    for(auto& s : coord_list) kernel.row(row_idx++) <<   0.     ,  0.     , (   s.x);   // edge 2 dof 5
                break;
            case 3: for(auto& s : coord_list) kernel.row(row_idx++) << -(1.-s.x), (1.-s.x),  0.;   // edge 3 dof 6
                    for(auto& s : coord_list) kernel.row(row_idx++) << -(   s.x), (   s.x),  0.;   // edge 3 dof 7
                break;
            case 4: for(auto& s : coord_list) kernel.row(row_idx++) << -(1.-s.x),  0., (1.-s.x);   // edge 4 dof 8
                    for(auto& s : coord_list) kernel.row(row_idx++) << -(   s.x),  0., (   s.x);   // edge 4 dof 9
                break;   
            case 5: for(auto& s : coord_list) kernel.row(row_idx++) <<   0.,-(1.-s.x), (1.-s.x);   // edge 5 dof 10
                    for(auto& s : coord_list) kernel.row(row_idx++) <<   0.,-(   s.x), (   s.x);   // edge 5 dof 11
                break;  
            }
        }else if(entity_dim == 2){  // face
            switch (entity_idx)
            {
            case 0: for(auto& s : coord_list) kernel.row(row_idx++) <<  1.       , 0.       , 0.       ;    // face 0 dof 12
                    for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 1.       , 0.       ;    // face 0 dof 13
                break;
            case 1: for(auto& s : coord_list) kernel.row(row_idx++) <<  1.       , 0.       , 0.       ;    // face 0 dof 14
                    for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 0.       , 1.       ;    // face 0 dof 15
                break;
            case 2: for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 1.       , 0.       ;    // face 0 dof 16
                    for(auto& s : coord_list) kernel.row(row_idx++) <<  0.       , 0.       , 1.       ;    // face 0 dof 17
                break;
            case 3: for(auto& s : coord_list) kernel.row(row_idx++) << -inv_sqrt3, inv_sqrt3, 0.       ;    // face 0 dof 18
                    for(auto& s : coord_list) kernel.row(row_idx++) << -inv_sqrt3, 0.       , inv_sqrt3;    // face 0 dof 19
                break;
            }
        }
        break;
    }
    default:
        Logger::warning("Hcurl_tetrahedron::get_dof_signature: higher order case not available.");
    }
}



void Hcurl_triangle::get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("Hcurl_triangle::get_basis_s - Basis functions in H(curl) are vector-valued, call get_basis_v instead.");
    throw std::logic_error("Basis functions in H(curl) are vector-valued, call get_basis_v instead.");
}

void Hcurl_triangle::get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("Hcurl_triangle::get_ED_basis_s - Exterior derivative of basis functions in H(curl) corresponds to the vector-valued curl, call get_ED_basis_v instead.");
    throw std::logic_error("Exterior derivative of basis functions in H(curl) corresponds to the vector-valued curl, call get_ED_basis_v instead.");
}

void Hcurl_tetrahedron::get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("Hcurl_tetrahedron::get_basis_s - Basis functions in H(curl) are vector-valued, call get_basis_v instead.");
    throw std::logic_error("Basis functions in H(curl) are vector-valued, call get_basis_v instead.");
}

void Hcurl_tetrahedron::get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("Hcurl_tetrahedron::get_ED_basis_s - Exterior derivative of basis functions in H(curl) corresponds to the vector-valued curl, call get_ED_basis_v instead.");
    throw std::logic_error("Exterior derivative of basis functions in H(curl) corresponds to the vector-valued curl, call get_ED_basis_v instead.");
}
