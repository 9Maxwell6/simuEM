#include "math/fem/space_H1.h"

#include <stdexcept>


using namespace simu;


H1_Space::H1_Space(int dim, int p) : FEM_Space(dim, p){}


bool H1_Space::add_basis_shape(Basis_Shape g)
{
    std::unique_ptr<H1_Space> shape_;
    switch (g) {
        case Basis_Shape::TETRAHEDRON: shape_ = std::unique_ptr<H1_Space>(new H1_tetrahedron(p_)); break;
        default: 
        {
            Logger::warning("H1_Space::add_geometry - geometry not available: return false");
            return false;
        }
    }

    // uniqueness check
    for (auto& existing : shape_H1_)
    {
        if (typeid(*existing.second) == typeid(*shape_)) return false;
    }

    //shape_H1_.push_back(std::move(shape_));
    shape_H1_[g] = std::move(shape_);
    basis_shapes_.push_back(g);
    return true;
}


FEM_Space * H1_Space::get_basis_space(Basis_Shape s) const
{
    auto it = shape_H1_.find(s);
    if (it != shape_H1_.end()) return it->second.get();

    Logger::error("H1_tetrahedron::get_basis_space - failed: key not found, return nullptr.");
    return nullptr;
}


const std::vector<Basis_Shape>& H1_Space::get_basis_shapes() const
{
    return basis_shapes_;
}



//void H1_Space::get_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const { shape_H1_.at(s)->get_basis_s(p, basis); }
//void H1_Space::get_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const { shape_H1_.at(s)->get_basis_v(p, basis); }
//void H1_Space::get_ED_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const { shape_H1_.at(s)->get_ED_basis_s(p, basis); }
//void H1_Space::get_ED_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const { shape_H1_.at(s)->get_ED_basis_v(p, basis); }


/**
 * ================================= Triangle =================================
 * 
 */

H1_triangle::H1_triangle(int p) : H1_Space(2, p)
{
    n_node_   = 3;
    n_edge_   = 3;
    n_face_   = 1;
    n_volume_ = 0;


    n_dof_            = (p+1)*(p+2)/2;
    n_dof_per_node_   = 1;
    n_dof_per_edge_   = (p-1);
    n_dof_per_face_   = (p-1)*(p-2)/2;
    n_dof_per_volume_ = 0;
}

void H1_triangle::get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    switch(p_)
    {
        case 1:
            //  Barycentric coordinates λ
            //          λ0 = 1.0 - x - y;
            //          λ1 = x;
            //          λ2 = y;
            //
            if (basis.size() != 3) {
                throw std::invalid_argument("vector must be 3x1 for p-1 H(grad) Triangle.");
            }

            basis(0) = 1.0 - x - y;      // Vertex 0
            basis(1) = x;                // Vertex 1
            basis(2) = y;                // Vertex 2
            break;
        default:
            throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(p_));
    }
};

void H1_triangle::get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    switch(p_)
    {
        case 1:
            // grad(λ0) = grad(1-x-y-z) = [-1, -1]
            // grad(λ1) = grad(x)       = [ 1,  0]
            // grad(λ2) = grad(y)       = [ 0,  1]
            //
            if (basis.rows() != 3 || basis.cols() != 2) {
                throw std::invalid_argument("matrix must be 3x2 for p-1 H(grad) Triangle.");
            }

            basis << -1.0, -1.0,   // Vertex 0
                      1.0,  0.0,   // Vertex 1
                      0.0,  1.0;   // Vertex 2
            break;
        default:
            throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(p_));
    }
};


/**
 * =============================== Tetrahedron ===============================
 * 
 */

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



void H1_tetrahedron::get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;
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


void H1_tetrahedron::get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    double x = coord.x;
    double y = coord.y;
    double z = coord.z;
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




void H1_triangle::get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    Logger::error("H1_triangle::get_basis_v - Basis functions in H(grad) are scalar-valued, call get_basis_s instead.");
    throw std::logic_error("Basis functions in H(grad) are scalar-valued, call get_basis_s instead.");
};

void H1_triangle::get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("H1_triangle::get_ED_basis_s - Exterior derivative of basis functions in H(grad) corresponds to the vector-valued grad, call get_ED_basis_v instead.");
    throw std::logic_error("Exterior derivative of basis functions in H(grad) corresponds to the vector-valued grad, call get_ED_basis_v instead.");
}

void H1_tetrahedron::get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {
    Logger::error("H1_tetrahedron::get_basis_v - Basis functions in H(grad) are scalar-valued, call get_basis_s instead.");
    throw std::logic_error("Basis functions in H(grad) are scalar-valued, call get_basis_s instead.");
};

void H1_tetrahedron::get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {
    Logger::error("H1_tetrahedron::get_ED_basis_s - Exterior derivative of basis functions in H(grad) corresponds to the vector-valued grad, call get_ED_basis_v instead.");
    throw std::logic_error("Exterior derivative of basis functions in H(grad) corresponds to the vector-valued grad, call get_ED_basis_v instead.");
};
