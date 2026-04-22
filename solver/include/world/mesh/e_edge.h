#pragma once

#include "element.h"
namespace simu {

class Mesh;


/*-------------------------------------------------------------
 * local node id: 0, 1
 *                      
 *                      node_1, node_2
 * local edge id: 0 -> [0,      1     ]
 * 
 */
class Edge : public Element
{
protected:
    size_t node_idx_[2];

    
public:
    Edge() = default;
    Edge(const size_t* node_idx, size_t id, size_t property_id=0, int o=1);
    Edge(std::vector<std::size_t> node_idx, size_t id, size_t property_id=0, int o=1);
    Edge(size_t idx_1, size_t idx_2, size_t id, size_t property_id=0, int o=1);


    ~Edge(){};

    Geometry get_geometry() const override;

    const size_t * get_node_idx() const override;
    int get_node_num() const override{ return 2; }

    int get_geometry_node_num() const override { return o_ + 1; };

    int get_element_dim() const override { return 1; }

    void set_node_idx(const size_t *idx);
    void set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3);
    
    bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const override;

    void compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const override;
    void compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const override;
    
    void compute_dof_transformation_H_curl(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const override;
    void compute_dof_transformation_H_div (const Mesh& mesh, Eigen::Ref<MatrixXd> P) const override;
};


}