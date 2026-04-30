#pragma once

#include "math/field/field.h"

#include "math/fem/space_collection.h"
//#include "math/fem/element_data.h"


namespace simu {



class S_Field_fespace : public Field
{

private:
    size_t dof_counter = 0;
    const FEM_Space& fe_space_;
    const std::vector<dof_idx>& dof_;

    const std::vector<scalar_t>& value_;
    

public:
    S_Field_fespace(const Mesh& mesh,
                    const FEM_Space& fe_space, 
                    const std::vector<dof_idx>& dof, 
                    const std::vector<scalar_t>& value);

    double eval(const Ref_Coord& ref_coord, const Element& e) const;

    Field_Type get_field_type() const override { return Field_Type::FEM_SPACE; }

    void reset_dof_counter() {dof_counter = 0;};


};



class V_Field_fespace : public Field
{

private:
    size_t dof_counter = 0;
    const FEM_Space& fe_space_;
    const std::vector<dof_idx>& dof_;

    const std::vector<scalar_t>& value_;
    

public:
    V_Field_fespace(const Mesh& mesh, 
                    const FEM_Space& fe_space, 
                    const std::vector<dof_idx>& dof, 
                    const std::vector<scalar_t>& value);

    void eval(const Ref_Coord& ref_coord, const Element& e, Eigen::Ref<VectorXd> value) const;

    Field_Type get_field_type() const override { return Field_Type::FEM_SPACE; }

    void reset_dof_counter() {dof_counter = 0;};


};



}