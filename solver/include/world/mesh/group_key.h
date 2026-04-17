#pragma once

namespace simu {

struct Key 
{
    uint32_t dim;
    uint32_t id; 

    bool operator==(const Key& k) const {
        return dim == k.dim && id == k.id;
    }

    struct Hash 
    {
        size_t operator()(const Key& k) const 
        {
            return (static_cast<size_t>(k.dim) << 32) | k.id;
        }
    };

};

}