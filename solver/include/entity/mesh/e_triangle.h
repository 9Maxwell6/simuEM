#pragma once

#include "element.h"

namespace simu {

class Mesh;


/*-------------------------------------------------------------
 * local node id: 0, 1, 2
 *                      
 *                      node_1, node_2
 * local edge id: 0 -> [0,      1     ]
 *                1 -> [0,      2     ]
 *                2 -> [1,      2     ]
 *
 *                      node_1, node_2, node_3
 * local face id: 0 -> [0,      1,      2     ]

 */

class Triangle : public Element
{
protected:
    size_t node_idx_[3];

    
public:
    Triangle() = default;
    Triangle(const size_t* node_idx, size_t id, size_t property_id=0, int o=1);
    Triangle(std::vector<std::size_t> node_idx, size_t id, size_t property_id=0, int o=1);
    Triangle(size_t idx_1, size_t idx_2, size_t idx_3, size_t id, size_t property_id=0, int o=1);


    ~Triangle(){};

    Geometry get_geometry() const override;

    const size_t * get_node_idx() const override;
    int get_node_num() const override{ return 3; }


    void set_node_idx(const size_t *idx);
    void set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3);
    
    bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const override;
    void compute_D_shape(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const override;

    void compute_dof_transformation_H_curl(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const override;
    void compute_dof_transformation_H_div (const Mesh& mesh, Eigen::Ref<MatrixXd> P) const override;
};




}