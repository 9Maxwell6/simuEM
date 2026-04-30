#include "math/field/field_fespace.h"


using namespace simu;

S_Field_fespace::S_Field_fespace(const Mesh& mesh, 
                                 const FEM_Space& fe_space, 
                                 const std::vector<dof_idx>& dof, 
                                 const std::vector<scalar_t>& value) : Field(mesh), fe_space_(fe_space), dof_(dof), value_(value)
{
    //TODO
}

V_Field_fespace::V_Field_fespace(const Mesh& mesh, 
                                 const FEM_Space& fe_space, 
                                 const std::vector<dof_idx>& dof, 
                                 const std::vector<scalar_t>& value) : Field(mesh), fe_space_(fe_space), dof_(dof), value_(value)
{
    //TODO
}

double S_Field_fespace::eval(const Ref_Coord& ref_coord, const Element& e) const
{
    Logger::error("S_Field_fespace::eval - not implemented.");
    return 0;
}


void V_Field_fespace::eval(const Ref_Coord& ref_coord, const Element& e, Eigen::Ref<VectorXd> value) const
{
    Logger::error("S_Field_fespace::eval - not implemented.");
}
