#pragma once

#include "math/operator/operator.h"
#include "math/fem/shape.h"

#include "math/data_format.h"






namespace simu {


class Interpolator : public Operator
{

protected:
    Interpolator();


public:
    static constexpr int SIZE = 1;

};


/**
 * interpolation from nodal(N) dof to edge(E) dof.
 */
class Interpolator_H1_to_Hcurl : public Interpolator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_1 && s_2 != Space::H_curl)
        {
            Logger::error("Interpolator_H1_to_Hcurl: require H_1 -> H_curl.");
            throw std::invalid_argument("Interpolator_H1_to_Hcurl: require H_1 -> H_curl.");
        }
    }

public:
    static constexpr int INTERPOLATOR_ID = 0;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
    
};


}