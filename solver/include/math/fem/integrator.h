#pragma once

#include "math/fem/integration.h"


namespace simu {





/**
 * s: scalar coefficient
 * v: vector coefficient
 * S: scalar function
 * V: vector function
 */
class Integrator
{

protected:
    Integrator();

    std::vector<Integration_Point *> Triangle_order_;
    std::vector<Integration_Point *> Tetrahedron_order_;
    // TODO： support more integration rules: point/edge/squares/cube/ etc.s

    virtual void assemble_element_matrix();
    

//public:

    
};



/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_S__S : public Integrator
{
    void assemble_element_matrix() override;
};


/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_grad_S__grad_S : public Integrator
{
    void assemble_element_matrix() override;
};


/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_V__grad_S : public Integrator
{
    void assemble_element_matrix() override;
};


/**
 * domain: H_curl
 * range:  H_curl
 * 
 */
class Integrator__s_curl_V__curl_V : public Integrator
{
    void assemble_element_matrix() override;
};

/**
 * domain: H_curl, H_div
 * range:  H_curl, H_div
 * 
 */
class Integrator__s_V__V : public Integrator
{
    void assemble_element_matrix() override;
};


}