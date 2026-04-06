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
 * All bilinear-form integrators.
 * 
 * 
 * 
 */


/**
 * domain: H_1
 * range:  H_1 
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s_S__S : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_1 || s_2 != Space::H_1)
        {
            Logger::error("Integrator__s_S__S: require H_1 -> H_1.");
            throw std::invalid_argument("Integrator__s_S__S: require H_1 -> H_1.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 0;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
    
};




/**
 * domain: H_1
 * range:  H_1 
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s_grad_S__grad_S : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_1 || s_2 != Space::H_1)
        {
            Logger::error("Integrator__s_grad_S__grad_S: require H1 -> H1.");
            throw std::invalid_argument("Integrator__s_grad_S__grad_S: require H1 -> H1.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 1;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


/**
 * domain: H_1
 * range:  H_1 
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator_H1__s_V__grad_S : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_1 || s_2 != Space::H_1)
        {
            Logger::error("Integrator_H1__s_V__grad_S: require H_1 -> H_1.");
            throw std::invalid_argument("Integrator_H1__s_V__grad_S: require H_1 -> H_1.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 2;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


/**
 * domain: H_curl
 * range:  H_curl
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s_curl_V__curl_V : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if(s_1 != Space::H_curl || s_2 != Space::H_curl)
        {
            Logger::error("Integrator__s_curl_V__curl_V: require H_curl -> H_curl.");
            throw std::invalid_argument("Integrator__s_curl_V__curl_V: require H_curl -> H_curl.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 3;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};

/**
 * domain: H_curl, H_div
 * range:  H_curl, H_div
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s_V__V : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if((s_1 != Space::H_curl || s_2 != Space::H_curl) && (s_1 != Space::H_div || s_2 != Space::H_div) )
        {
            Logger::error("Integrator__s_V__V: require H_curl -> H_curl or H_div -> H_div.");
            throw std::invalid_argument("Integrator__s_V__V: require H_curl -> H_curl or H_div -> H_div.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 4;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


/**
 * domain: H_curl, H_div
 * range:  H_1
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 * TODO: currently only support Hcurl-H1, extends to Hdiv-H1.
 */
class Integrator__s_V__grad_S : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if((s_1 != Space::H_curl || s_2 != Space::H_1) && (s_1 != Space::H_div || s_2 != Space::H_1) )
        {
            Logger::error("Integrator__s_V__grad_S: require H_curl -> H_1 or H_div -> H_1.");
            throw std::invalid_argument("Integrator__s_V__grad_S: require H_curl -> H_1 or H_div -> H_1.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 5;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};

/**
 * domain: H_1
 * range:  H_curl, H_div
 * 
 * @param coeff scalar coefficient value at each element.
 * @param e_data element data struct contains all information needed to assemble element matrix.
 * @param element_matrix computed local element matrix will be add to element_matrix.
 * 
 */
class Integrator__s_grad_S__V : public Integrator
{
private:
    void static check_precondition(const FEM_Space* space_1, const FEM_Space* space_2)
    {
        Space s_1 = space_1->get_function_space();
        Space s_2 = space_2->get_function_space();
        if((s_1 != Space::H_1 || s_2 != Space::H_curl) && (s_1 != Space::H_1 || s_2 != Space::H_div) )
        {
            Logger::error("Integrator__s_grad_S__V: require H_1 -> H_curl or H_1 -> H_div.");
            throw std::invalid_argument("Integrator__s_grad_S__V: require H_1 -> H_curl or H_1 -> H_div.");
        }
    }

public:
    static constexpr int INTEGRATOR_ID = 6;

    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};





}