#pragma once
#include "../property/property.h"
#include "math/data_format.h"
#include "math/ref_coord.h"
#include "utils/logger.h"
#include "entity/mesh/node.h"



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

    inline size_t get_id() const {return id_;}
    inline int    get_geometry_order() const {return o_;}
    inline size_t get_property_id() const {return property_id_;}
    inline void   set_property_id(size_t property_id) {property_id_ = property_id;}
    

    /**
     * @brief compute Jacobian matrix of element at each integration point (transformation from reference element to actual element in mesh).
     * 
     * @param mesh mesh that contains the element.
     * @param i_p integration point from quadrature rules.
     * @param J Jacobian matrix to be filled.
     * @return true if the Jacobian matrix is independent from integration points, this is for avoiding recompute same Jacobian matrix.
     */
    virtual bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const = 0;

    virtual void compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const = 0;
    virtual void compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const = 0;

    /**
     * @brief compute transformation matrix to correct dof direction on reference element to actual element in H_curl space.
     * 
     * @param mesh mesh that contains the element.
     * @param P transformation matrix to be filled.
     */
    virtual void compute_dof_transformation_H_curl(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const = 0;

    /**
     * @brief compute transformation matrix to correct dof direction on reference element to actual element in H_div space.
     * 
     * @param mesh mesh that contains the element.
     * @param P transformation matrix to be filled.
     */
    virtual void compute_dof_transformation_H_div(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const = 0;

};

}