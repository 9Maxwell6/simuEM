#pragma once

#include "element.h"

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

    Type get_Type() const override;

    const size_t * get_NodeIdx() const override;


    void set_NodeIdx(const size_t *idx);
    void set_NodeIdx(size_t idx_1, size_t idx_2, size_t idx_3);
    
    
    

};