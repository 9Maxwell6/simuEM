#pragma once

#include "fem_space.h"


namespace simu {



class Hcurl_Space : public FEM_Space 
{
private:
    std::vector<std::unique_ptr<Hcurl_Space>> shape_Hcurl_;
    
public:
    Hcurl_Space(int dim, int p = 1);

    bool add_basis_shape(Basis_Shape g) override;

    // Returns the number of DOFs per element
    virtual int get_element_dof() const = 0;

    // Returns basis values at a point in the unit tetrahedron
    // For H1: Scalars. For HCurl: Vectors.
    virtual void get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const = 0;
    virtual void get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const = 0;

    // Return vector proxy of Exterior Derivative of basis of the corresponding form.
    virtual void get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const = 0;
    virtual void get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const = 0;

    Space get_function_space() const override {return Space::H_curl;};
};


/**
 * Reference tetrahedron:
 * 
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0]  a3=[0,0,1]    
 * 
 */
class Hcurl_tetrahedron : public Hcurl_Space {
    friend class Hcurl_Space;
private:
    std::vector<std::vector<double>> basis_cache;

    Hcurl_tetrahedron(int p = 1);
    ~Hcurl_tetrahedron(){};

public:
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