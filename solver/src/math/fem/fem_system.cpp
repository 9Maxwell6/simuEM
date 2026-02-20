#include "math/fem/fem_system.h"

using namespace simu;



FEM_System::FEM_System(Mesh& mesh):mesh_(mesh){};



/**
 * @brief Assign space to all elements in mesh
 *
 * @param fs finite element space.
 * @param p_order order of basis function.
 */
void FEM_System::assign_FE_space(Space fs, int p_order)
{
    // TODO: should use vector<Space>, what if we assign duplicate Space?
    // use vector<FEM_Space *>,
    //global_space.push_back(fs);
}
