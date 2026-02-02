#pragma once

#include <stddef.h>
#include <vector>
#include "info.h"


enum Type { 
    POINT, 
    SEGMENT, 
    TRIANGLE, 
    TETRAHEDRON
    // currently not support other shape
};







class Element
{
protected:
    size_t property_id_;    // id for element info

public:

    virtual ~Element() {}

    /// Returns element's type
    virtual Type get_Type() const = 0;

    virtual const size_t * get_NodeIdx() const= 0;
    virtual void set_NodeIdx(const size_t *ind) = 0;

    inline size_t get_propertyId() const {return property_id_;}
    inline void   set_propertyId(size_t property_id) {property_id_ = property_id;}



};