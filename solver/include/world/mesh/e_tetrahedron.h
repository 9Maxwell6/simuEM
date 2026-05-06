#pragma once
#include "element.h"

namespace simu {

class Mesh;


/*-------------------------------------------------------------
 * local node id: 0, 1, 2, 3
 *                      
 *                      node_1, node_2
 * local edge id: 0 -> [0,      1     ]
 *                1 -> [0,      2     ]
 *                2 -> [0,      3     ]
 *                3 -> [1,      2     ]
 *                4 -> [1,      3     ]
 *                5 -> [2,      3     ]
 * 
 *                      node_1, node_2, node_3
 * local face id: 0 -> [0,      1,      2     ]
 *                1 -> [0,      1,      3     ]
 *                2 -> [0,      2,      3     ]
 *                3 -> [1,      2,      3     ]
 * 
 * 
 * Reference tetrahedron:
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0]  a3=[0,0,1]    
 *
 */

class Tetrahedron : public Element
{
protected:
    size_t  node_idx_[4];

    
public:
    Tetrahedron() = default;
    Tetrahedron(const size_t* node_idx, size_t id, size_t property_id=0, int o=1);
    Tetrahedron(std::vector<std::size_t> node_idx, size_t id, size_t property_id=0, int o=1);
    Tetrahedron(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4, size_t id, size_t property_id=0, int o=1);


    ~Tetrahedron(){};

    Geometry get_geometry() const override;

    const size_t * get_node_idx() const override;
    int get_node_num() const override{ return 4; }

    int get_geometry_node_num() const override { return (o_ + 1) * (o_ + 2) * (o_ + 3) / 6; };

    int get_dim() const override { return 3; }
    int get_n_node() const override { return 4; }
    int get_n_edge() const override { return 6; }
    int get_n_face() const override { return 4; }
    int get_n_cell() const override { return 1; }

    void set_node_idx(const size_t *idx);
    void set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4);
    
    bool compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const override;



    void compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const override;
    void compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const override;

    std::vector<Ref_Coord> edge_map(const std::vector<Integration_Point>& edge_coord, size_t edge_idx) const override;
    std::vector<Ref_Coord> face_map(const std::vector<Integration_Point>& face_coord, size_t face_idx) const override;

    void tangent(Eigen::Ref<MatrixXd> t) const override;
    void normal(Eigen::Ref<MatrixXd> n) const override;

};

}