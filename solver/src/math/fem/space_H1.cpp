#include "math/fem/space_H1.h"

#include <stdexcept>


using namespace simu;


H1_Space::H1_Space(int dim, int p) : FEM_Space(dim, p){}


H1_tetrahedron::H1_tetrahedron(int p) : H1_Space(3, p)
{
    n_node_   = 4;
    n_edge_   = 6;
    n_face_   = 4;
    n_volume_ = 1;


    n_dof_            = (p+1)*(p+2)*(p+3)/6;
    n_dof_per_node_   = 1;
    n_dof_per_edge_   = (p-1);
    n_dof_per_face_   = (p-1)*(p-2)/2;
    n_dof_per_volume_ = (p-1)*(p-2)*(p-3)/6;
}



void H1_tetrahedron::get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const {
    double x = p.x;
    double y = p.y;
    double z = p.z;
    switch(p_)
    {
        case 1:
            //  Barycentric coordinates λ
            //          λ0 = 1.0 - x - y - z;
            //          λ1 = x;
            //          λ2 = y;
            //          λ3 = z;
            //
            if (basis.size() != 4) {
                throw std::invalid_argument("vector must be 4x1 for p-1 H(grad) Tetrahedron.");
            }

            basis(0) = 1.0 - x - y - z;  // Vertex 0
            basis(1) = x;                // Vertex 1
            basis(2) = y;                // Vertex 2
            basis(3) = z;                // Vertex 3
            break;
        default:
            throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(p_));
    }
};


void H1_tetrahedron::get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const {
    double x = p.x;
    double y = p.y;
    double z = p.z;
    switch(p_)
    {
        case 1:
            // grad(λ0) = grad(1-x-y-z) = [-1, -1, -1]
            // grad(λ1) = grad(x)       = [ 1,  0,  0]
            // grad(λ2) = grad(y)       = [ 0,  1,  0]
            // grad(λ3) = grad(z)       = [ 0,  0,  1]
            //
            if (basis.rows() != 4 || basis.cols() != 3) {
                throw std::invalid_argument("matrix must be 4x3 for p-1 H(grad) Tetrahedron.");
            }

            basis << -1.0, -1.0, -1.0,  // Vertex 0
                      1.0,  0.0,  0.0,  // Vertex 1
                      0.0,  1.0,  0.0,  // Vertex 2
                      0.0,  0.0,  1.0;  // Vertex 3
            break;
        default:
            throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(p_));
    }

    

};





int H1_tetrahedron::get_element_dof() const{
    switch(p_)
    {
        case 1: return 4;
        default:
            throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(p_));
    }
}



void H1_tetrahedron::get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const {
    throw std::logic_error("Basis functions in H(grad) are scalar-valued, call get_basis_s instead.");
};

void H1_tetrahedron::get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const {
    throw std::logic_error("Exterior derivative of basis functions in H(grad) corresponds to the vector-valued grad, call get_ED_basis_v instead.");
};
