#pragma once

namespace simu {


/**
 * s: scalar coefficient
 * v: vector coefficient
 * S: scalar function
 * V: vector function
 */
class Integrator
{


    
};



/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_S__S : public Integrator
{

};


/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_grad_S__grad_S : public Integrator
{

};


/**
 * domain: H_1
 * range:  H_1 
 * 
 */
class Integrator__s_V__grad_S : public Integrator
{

};


/**
 * domain: H_curl
 * range:  H_curl
 * 
 */
class Integrator__s_curl_V__curl_V : public Integrator
{

};

/**
 * domain: H_curl, H_div
 * range:  H_curl, H_div
 * 
 */
class Integrator__s_V__V : public Integrator
{

};


}