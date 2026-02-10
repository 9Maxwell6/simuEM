//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"

#include <Eigen/Dense>

#include "math/fem/integration.h"
#include "math/fem/fem_space.h"
#include "math/fem/space_H1.h"
#include "math/fem/space_Hcurl.h"


#include <stdio.h>
//#include <config.h>
#include "entity/mesh/mesh_parser.h"


int main() {
    Mesh_Parser mp(Mesh_Format::GMSH);
    Mesh mesh = mp.load_mesh(SCRIPT_PATH "test_mesh_0.msh");
    return 0;
}