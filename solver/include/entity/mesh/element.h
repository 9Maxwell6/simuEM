#pragma once
#include "../property/property.h"


#include <stddef.h>
#include <vector>


enum class Type { 
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
    virtual Type get_Type() const = 0;

    virtual const size_t * get_nodeIdx() const = 0;
    virtual int get_nodeNum() const = 0;
    virtual void set_nodeIdx(const size_t *ind) = 0;

    inline size_t get_Id() const {return id_;}
    inline size_t get_propertyId() const {return property_id_;}
    inline void   set_propertyId(size_t property_id) {property_id_ = property_id;}



};