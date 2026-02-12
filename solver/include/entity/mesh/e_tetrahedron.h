#pragma once

#include "element.h"

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
 * local face id: 0 -> [1,      3,      2     ]
 *                1 -> [2,      3,      0     ]
 *                2 -> [1,      3,      0     ]
 *                3 -> [0,      1,      2     ]
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

    Type get_Type() const override;

    const size_t * get_nodeIdx() const override;
    int get_nodeNum() const override{ return 4; }


    void set_nodeIdx(const size_t *idx);
    void set_nodeIdx(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4);
    
    
    

};