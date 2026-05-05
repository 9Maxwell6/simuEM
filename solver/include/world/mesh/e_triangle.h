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

    int get_geometry_node_num() const override { return (o_ + 1) * (o_ + 2) / 2; };

    int get_dim() const override { return 2; }
    int get_n_node() const override { return 3; }
    int get_n_edge() const override { return 3; }
    int get_n_face() const override { return 1; }
    int get_n_volume() const override { return 0; }

    void set_node_idx(const size_t *idx);
    void set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3);
    
    bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const override;

    void compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const override;
    void compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const override;

    std::vector<Ref_Coord> edge_map(const std::vector<Integration_Point>& edge_coord, size_t edge_idx) const override;
    std::vector<Ref_Coord> face_map(const std::vector<Integration_Point>& face_coord, size_t face_idx) const override;

    void tangent(Eigen::Ref<MatrixXd> t) const override;
    void normal(Eigen::Ref<MatrixXd> n) const override;
};




}