#pragma once

#include "math/fem/operation.h"
#include "math/fem/integration.h"
#include "math/fem/integrator.h"
#include "math/data_format.h"
#include "math/fem/fem_util.h"



#include <Eigen/Dense>



namespace simu {


template<int phy_dim, int ref_dim>
struct Element_Data;


/**
 * All linear-form integrators.
 * 
 * 
 * 
 */



/**
 * space: H_1
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s__S : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        if(space_1 != nullptr && space_2 == nullptr)
        {
            Space s_1 = space_1->get_function_space();
            if(s_1 != Space::H_1)
            {
                Logger::error("Integrator__s__S: require H_1.");
                throw std::invalid_argument("Integrator__s__S: require H_1.");
            }
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 0;

    template<typename Field, int phy_dim, int ref_dim, typename Vec_Type>
    void static assemble_element_vector(Field&& F, Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_matrix);
    
};



}