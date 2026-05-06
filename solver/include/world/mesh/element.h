#pragma once
#include "../property/property.h"
#include "math/data_format.h"
#include "math/ref_coord.h"
#include "utils/logger.h"
#include "world/mesh/node.h"



#include <stddef.h>
#include <vector>

namespace simu {


enum class Geometry { 
    NODE, 
    EDGE, 
    TRIANGLE, 
    TETRAHEDRON
    // currently not support other shape
};





class Mesh;

class Element
{
protected:
    size_t id_;             // element id.
    size_t property_id_;    // id for element info.
    int o_;                 // geometry order of the actual element.

public:
    Element() = default;
    Element(size_t id, size_t property_id, int o): id_(id), property_id_(property_id), o_(o){};
    virtual ~Element() {}

    /// Returns element's type
    virtual Geometry get_geometry() const = 0;

    virtual const size_t * get_node_idx() const = 0;
    virtual int get_node_num() const = 0; // number of vertices
    virtual void set_node_idx(const size_t *ind) = 0;

    virtual int get_geometry_node_num() const = 0; // number of node that define the geometry.

    virtual int get_dim() const = 0;  // dimension of reference element.
    virtual int get_n_node() const = 0;
    virtual int get_n_edge() const = 0;
    virtual int get_n_face() const = 0;
    virtual int get_n_cell() const = 0;

    inline size_t get_id() const {return id_;}
    inline int    get_geometry_order() const {return o_;}
    inline size_t get_property_id() const {return property_id_;}
    inline void   set_property_id(size_t property_id) {property_id_ = property_id;}


    /**
     * @brief compute the physical coordinate of the reference coordinate of the element. 
     * 
     * @param mesh mesh that contains the element.
     * @param ref_coord reference coordinate in the element.
     * @param phy_coord physical coordinate in mesh to be computed.
     */
    void physical_point(const Mesh& mesh, const Ref_Coord& ref_coord, Eigen::Ref<VectorXd> phy_coord) const;
    

    /**
     * @brief compute Jacobian matrix of element at each integration point (transformation from reference element to actual element in mesh).
     * 
     * @param mesh mesh that contains the element.
     * @param i_p integration point from quadrature rules.
     * @param J Jacobian matrix to be filled.
     * @return true if the Jacobian matrix is independent from integration points, this is for avoiding recompute same Jacobian matrix.
     */
    virtual bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const = 0;


    /**
     * @brief Reference-space tangent vectors of a face.
     *
     * usually combined with the volume Jacobian J at a face quadrature point, 
     * they give the face Jacobian columns J*e1, J*e2 and hence
     * the oriented area vector
     *
     *     a_f = (J e1) x (J e2).
     *
     *
     * @param[in]  face_idx local face index, in `[0, num_faces())`.
     * @param[out] e1, e2   refeference tangent vectors; caller must pre-size to the
     *                      element's reference dimension.
     */
    virtual void face_ref_edge(int face_idx, Eigen::Ref<VectorXd> e1, Eigen::Ref<VectorXd> e2) const = 0;

    virtual void compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const = 0;
    virtual void compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const = 0;

    /**
     * @brief convert reference coordinate defined on edge to the coordinate of each edge of the reference element.
     * 
     * @param edge_coord list of reference coordinate on edge.
     * @return vector of reference coordinate on one edge of the reference element.
     */
    virtual std::vector<Ref_Coord> edge_map(const std::vector<Integration_Point>& edge_coord, size_t edge_idx) const = 0;

    /**
     * @brief convert reference coordinate defined on face to the coordinate of each face of the reference element.
     * 
     * @param face_coord list of reference coordinate on edge.
     * @return vector of reference coordinate on one face of the reference element, 
     */
    virtual std::vector<Ref_Coord> face_map(const std::vector<Integration_Point>& face_coord, size_t face_idx) const = 0;

    /**
     * @brief compute the tangent vector along each edge, with magnitude equal to the length of the reference edge,
     *        and store them to output matrix.
     * 
     * @param t output matrix, with each row represnets one tangent on one edge of the reference element, 
     */
    virtual void tangent(Eigen::Ref<MatrixXd> t) const = 0;


    /**
     * @brief compute the outward normal vector of each surface, with magnitude equal to the area of the reference face,
     *        and store them to output matrix.
     * 
     * @param n output matrix, with each row represnets one tangent on one edge of the reference element, 
     * 
     * @note for 2D elements, like triangle, normal is normal vector of each edge.
     */
    virtual void normal(Eigen::Ref<MatrixXd> n) const = 0;
    
    

};

}