#pragma once
#include "../property/property.h"


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
    virtual int get_node_num() const = 0;
    virtual void set_node_idx(const size_t *ind) = 0;

    inline size_t get_id() const {return id_;}
    inline int    get_geometry_order() const {return o_;}
    inline size_t get_property_id() const {return property_id_;}
    inline void   set_property_id(size_t property_id) {property_id_ = property_id;}
    



};

}