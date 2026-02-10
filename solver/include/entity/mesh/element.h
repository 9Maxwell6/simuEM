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
    size_t property_id_;    // id for element info
    int o_;                 // geometry order of the actual element.

public:
    Element(size_t property_id, int o): property_id_(property_id), o_(o){};
    virtual ~Element() {}

    /// Returns element's type
    virtual Type get_Type() const = 0;

    virtual const size_t * get_NodeIdx() const= 0;
    virtual void set_NodeIdx(const size_t *ind) = 0;

    inline size_t get_propertyId() const {return property_id_;}
    inline void   set_propertyId(size_t property_id) {property_id_ = property_id;}



};