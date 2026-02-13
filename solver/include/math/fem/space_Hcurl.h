#pragma once

#include "fem_space.h"


namespace simu {


/**
 * Reference tetrahedron:
 * 
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0]  a3=[0,0,1]    
 * 
 */
class Hcurl_tetrahedron : public FEM_Space {
private:
    std::vector<std::vector<double>> basis_cache;

public:
    Hcurl_tetrahedron(int p = 1);

    ~Hcurl_tetrahedron(){};

    //void get_basis_v(Integration_Point p, const Eigen::Ref<MatrixXd> basis) const override {
    //    const double * ss = basis.data();
    //};


    void get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const override;
    void get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const override;

    void get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const override;
    void get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const override;

    int get_element_dof() const override;
    Space get_function_space() const override {return Space::H_curl;};

};


}