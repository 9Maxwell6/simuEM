#pragma once

#include "entity/mesh/mesh.h"
#include "math/fem/integration.h"
#include "math/data_format.h"
#include "utils/logger.h"


namespace simu {


template <typename Ref_Element, int phy_dim>
struct Element_Transformation
{
    const Mesh* mesh_ = nullptr;
    static constexpr int ref_dim  = Ref_Element::ref_dim;
    static constexpr int node_num = Ref_Element::node_num;

    Matrix<phy_dim, ref_dim>  J_;
    Matrix<ref_dim, phy_dim>  inv_J_;
    double det_J_;

    Element_Transformation() = default;
    void set_mesh(const Mesh& mesh) { mesh_ = &mesh; }

    Element_Transformation(const Mesh& mesh): mesh_(&mesh){};



    const Matrix<phy_dim, Ref_Element::ref_dim>& compute_Jacobian(const Integration_Point& i_p, const Element * e, const Matrix<phy_dim, node_num>& node_matrix);

    const Matrix<phy_dim, Ref_Element::ref_dim>& compute_Jacobian(const Element * e);
    
    const Matrix<Ref_Element::ref_dim, phy_dim>& compute_inv_Jacobian();

    double compute_det_Jacobian();
    

};


struct Triangle_o1 
{
    static constexpr int ref_dim  = 2;
    static constexpr int node_num = 3;
    static constexpr bool constant_jacobian = true;

    static void D_shape(const Element * e, const Integration_Point& i_p, Matrix<3, 2>& d_shape);

    template <int phy_dim>
    static void Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, 2>& J);
};


struct Tetrahedron_o1 {
    static constexpr int ref_dim  = 3;
    static constexpr int node_num = 4;
    static constexpr bool constant_jacobian = true;

    static void D_shape(const Element * e, const Integration_Point& i_p, Matrix<4, 3>& d_shape);

    template <int phy_dim>
    static void Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, 3>& J);
};


struct General_Element {
    static constexpr int ref_dim  = Eigen::Dynamic;
    static constexpr int node_num = Eigen::Dynamic;
    static constexpr bool constant_jacobian = false;

    int runtime_ref_dim;
    int runtime_node_num;

    static void D_shape(const Element * e, const Integration_Point& i_p, MatrixXd& d_shape) {
        Geometry g = e->get_geometry();
        int o = e->get_geometry_order();
        // TODO: handle each geometry and shape order
    }

    template <int phy_dim>
    static void Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, ref_dim>& J);
};



using Transform_Triangle_o1_2D       = Element_Transformation<Triangle_o1, 2>;
using Transform_Triangle_o1_3D       = Element_Transformation<Triangle_o1, 3>;
using Transform_Tetrahedron_o1_3D    = Element_Transformation<Tetrahedron_o1, 3>;

using Transform_general_2D           = Element_Transformation<General_Element, 2>;
using Transform_general_3D           = Element_Transformation<General_Element, 3>;



}