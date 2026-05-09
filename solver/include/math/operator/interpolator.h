#pragma once

#include "math/operator/operator.h"
#include "math/fem/shape.h"

#include "math/data_format.h"

#include "math/dof/dof_manager.h"





namespace simu {


class Interpolator : public Operator
{

protected:
    Interpolator();


public:
    static constexpr int SIZE = 2;

};


/**
 * interpolation from nodal(N) dof to edge(E) dof.
 */
class Interpolator__H1_to_Hcurl : public Interpolator
{
private:
    void static do_once(const FEM_Space* space_1, const FEM_Space* space_2, DoF_Manager* dm, const std::array<size_t,4>* e_size=nullptr)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_curl && s_2 != Space::H_1)
            Logger::error("Interpolator__H1_to_Hcurl: require H_1(column) -> H_curl(row).");

        if(space_2->get_vdim() != space_2->get_dim())
            Logger::error("Interpolator__H1_to_Hcurl: require vdim = dim in H1 space.");
        
        if(dm){
            if(!dm->is_ready())
                if(e_size){
                    dm->initialize((*e_size)[0], (*e_size)[1], (*e_size)[2], (*e_size)[3]);
                }else{
                    dm->initialize();
                }
            
        }else{
            Logger::error("Interpolator__H1_to_Hcurl: missing hash table.");
        }
    }

public:
    static constexpr int INTERPOLATOR_ID = 0;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
    
};










/**
 * interpolation from nodal(N) dof to edge(E) dof.
 */
class Interpolator__grad_H1_to_Hcurl : public Interpolator
{
private:
    void static do_once(const FEM_Space* space_1, const FEM_Space* space_2, DoF_Manager* dm, const std::array<size_t,4>* e_size=nullptr)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_curl && s_2 != Space::H_1)
            Logger::error("Interpolator__H1_to_Hcurl: require H_1(column) -> H_curl(row).");

        if(space_2->get_vdim() != 1)
            Logger::error("Interpolator__H1_to_Hcurl: require scalar H1 space, vdim=1.");
        
        if(dm){
            if(!dm->is_ready())
                if(e_size){
                    dm->initialize((*e_size)[0], (*e_size)[1], (*e_size)[2], (*e_size)[3]);
                }else{
                    dm->initialize();
                }
            
        }else{
            Logger::error("Interpolator__H1_to_Hcurl: missing hash table.");
        }
    }

public:
    static constexpr int INTERPOLATOR_ID = 1;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
    
};






}