#pragma once

#include "element.h"

/*-------------------------------------------------------------
 * local node id: 0, 1, 2, 3
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
    size_t  node_idx_[3];

    
public:

    Triangle(const size_t* node_idx, size_t property_id=0, int o=1);
    Triangle(size_t idx_1, size_t idx_2, size_t idx_3, size_t property_id=0, int o=1);


    ~Triangle(){};

    Type get_Type() const override;

    const size_t * get_NodeIdx() const override;


    void set_NodeIdx(const size_t *idx);
    void set_NodeIdx(size_t idx_1, size_t idx_2, size_t idx_3);
    
    
    

};