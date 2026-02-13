#pragma once

#include "entity/mesh/mesh.h"
#include "math/fem/fem_space.h"
#include "utils/string_utils.h"


#include <functional>

namespace simu {


class T_Omega
{

private:
    Mesh& mesh_;

    int dim_;

    std::vector<Key> key_true_boundary;              // key to 1D/2D true boundary element groups
    std::vector<Key> key_conductor_interface;        // key to 1D/2D conductor boundary element groups
    std::vector<Key> key_source;                     // key to 2D/3D source element groups
    std::vector<Key> key_insulator;                  // key to 2D/3D insulating region element groups
    std::vector<Key> key_conductor;                  // key to 2D/3D conducting region element groups
    std::vector<Key> key_conductor_interface_layer;  // key to 2D/3D conducting region element groups in contact with conductor boundarys.
    


    /**
     * @brief Filter function used for Mesh::mark_elements, it will checks if 
     * at least one element in the conductor_interface group is covered by the target element
     *
     * element A is considered "covered" by element B if all A's vertices is contained in B's vertices
     *
     * @param e Pointer to the target element.
     * @return true  If at least one element in conductor_interface is covered by the target element.
     * @return false Otherwise
     */
    std::function<bool(Element*)> conductor_interface_layer_filter() {
        return [this](Element* e) -> bool {
            const size_t * e_ids = e->get_nodeIdx();
            int e_size = e->get_nodeNum();
            //std::cout<<e_ids[0]<<" "<< e_ids[1]<<" " << e_ids[2]<<" "<< e_ids[3]<<" " <<std::endl;
            for(Key& key: this->key_conductor_interface)
            {    
                for (Element* ie : this->mesh_.get_group(key))
                {
                    const size_t * ie_ids = ie->get_nodeIdx();
                    int ie_size = ie->get_nodeNum();

                    for (size_t i = 0; i < ie_size; ++i) 
                    {
                        for (size_t j = 0; j < e_size; ++j) 
                        {
                            if (ie_ids[i] == e_ids[j]) { return true; }
                        }
                    }
                    /*
                    bool match = true;
                    for (size_t i = 0; i < ie_size; ++i) 
                    {
                        bool found = false;
                        for (size_t j = 0; j < e_size; ++j) 
                        {
                            if (ie_ids[i] == e_ids[j]) { found = true; break; }
                        }
                        if (!found) { match = false; break; }
                    }
                    if (match) return true;
                    */
                }
            }
            return false;
        };
    }



public:
    T_Omega(Mesh& mesh);

};


}