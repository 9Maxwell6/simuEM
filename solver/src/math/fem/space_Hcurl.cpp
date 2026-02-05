#include "math/fem/space_Hcurl.h"

#include <stdexcept>

Hcurl_tetrahedron::Hcurl_tetrahedron(int p) : FEM_Space(p)
{
    dim_ = 3;

}






void Hcurl_tetrahedron::get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const 
{
    double x = p.x;
    double y = p.y;
    double z = p.z;

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

};

void Hcurl_tetrahedron::get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const {
    double x = p.x;
    double y = p.y;
    double z = p.z;
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
    
};






int Hcurl_tetrahedron::get_element_dof() const{
    switch(p_)
    {
        case 1: return 6;
        default:
            throw std::invalid_argument("Edge element not available for order:  "+std::to_string(p_));
    }
}


void Hcurl_tetrahedron::get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const {
    throw std::logic_error("Basis functions in H(curl) are vector-valued, call get_basis_v instead.");
};


void Hcurl_tetrahedron::get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const {
    throw std::logic_error("Exterior derivative of basis functions in H(curl) corresponds to the vector-valued curl, call get_ED_basis_v instead.");
};
