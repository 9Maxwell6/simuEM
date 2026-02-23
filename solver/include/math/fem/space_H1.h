#pragma once

#include "fem_space.h"


namespace simu {


class H1_Space : public FEM_Space 
{
private:
    std::vector<std::unique_ptr<H1_Space>> shape_H1_;

public:
    H1_Space(int dim, int p = 1);

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

    Space get_function_space() const override {return Space::H_1;};
};



/**
 * Reference tetrahedron:
 * 
 *      [x,y,z]
 *   a0=[0,0,0]  a1=[1,0,0]  a2=[0,1,0]  a3=[0,0,1]    
 *
 */
class H1_tetrahedron : public H1_Space 
{
    friend class H1_Space;
private:
    std::vector<std::vector<double>> basis_cache;

    H1_tetrahedron(int p = 1);
    ~H1_tetrahedron(){};

public:
    //void get_basis_s(Integration_Point p, const Eigen::Ref<VectorXd> basis) const override {
    //   const double * ss = basis.data();
    //};

    void get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const override;
    void get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const override;

    void get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const override;
    void get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const override;

    int get_element_dof() const override;
    
};


}