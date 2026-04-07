#pragma once

#include "math/field/field.h"

#include "math/fem/space_collection.h"

namespace simu {



template<int phy_dim>
class S_Field_fespace : public Field<phy_dim>
{

private:
    size_t dof_counter = 0;
    const FEM_Space& fe_space_;
    const std::vector<Element*>& elements_;
    const std::vector<dof_idx>& dof_;

    const std::vector<double>& value_;
    

public:
    S_Field_fespace(const FEM_Space& fe_space, 
                  const std::vector<Element*>& element, 
                  const std::vector<dof_idx>& dof, 
                  const std::vector<double>& value);

    template<int ref_dim>
    double eval(const Integration_Point& i_p, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //
    }

    void reset_dof_counter() {dof_counter = 0;};


};



template<int phy_dim>
class V_Field_fespace : public Field<phy_dim>
{

private:
    size_t dof_counter = 0;
    const FEM_Space& fe_space_;
    const std::vector<Element*>& elements_;
    const std::vector<dof_idx>& dof_;

    const std::vector<double>& value_;
    

public:
    V_Field_fespace(const FEM_Space& fe_space, 
                  const std::vector<Element*>& element, 
                  const std::vector<dof_idx>& dof, 
                  const std::vector<double>& value);

    template<int ref_dim>
    Vector<phy_dim> eval(const Integration_Point& i_p, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //
    }

    void reset_dof_counter() {dof_counter = 0;};


};



}