#pragma once

#include "math/field/field.h"

#include "math/fem/space_collection.h"

namespace simu {



class S_Field_fespace : public Field
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

    template<int phy_dim, int ref_dim>
    double eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //
    }

    void reset_dof_counter() {dof_counter = 0;};


};



class V_Field_fespace : public Field
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

    template<int phy_dim, int ref_dim>
    Vector<phy_dim> eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //
    }

    void reset_dof_counter() {dof_counter = 0;};


};



}