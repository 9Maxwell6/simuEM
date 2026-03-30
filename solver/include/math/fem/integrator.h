#pragma once

#include "math/fem/operation.h"
#include "math/fem/integration.h"
#include "math/matrix.h"
#include "math/fem/assemble_data.h"
#include "math/fem/fem_util.h"



#include <Eigen/Dense>



namespace simu {





/**
 * s: scalar coefficient
 * v: vector coefficient
 * S: scalar function
 * V: vector function
 */
class Integrator : public Operation
{

protected:
    Integrator();

    std::vector<Integration_Point *> Triangle_order_;
    std::vector<Integration_Point *> Tetrahedron_order_;
    // TODO： support more integration rules: point/edge/squares/cube/ etc.s

    bool require_J = false;
    bool require_inv_J = false;
    bool require_det_J = false;

    //virtual  void assemble_element_matrix();
public:
     


//public:

    
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
class Integrator__s_S__S : public Integrator
{
public:
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
public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_V__grad_S : public Integrator
{
public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


/**
 * domain: H_curl
 * range:  H_curl
 * 
 */
class Integrator__s_curl_V__curl_V : public Integrator
{

public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};

/**
 * domain: H_curl, H_div
 * range:  H_curl, H_div
 * 
 */
class Integrator__s_V__V : public Integrator
{
public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


}