//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"


#include "math/fem/integration.h"
#include "math/fem/fem_space.h"
#include "math/fem/space_H1.h"
#include "math/fem/space_Hcurl.h"

#include "utils/logger.h"

#include <Eigen/Dense>
#include <stdio.h>
//#include <config.h>
#include "entity/mesh/mesh_parser.h"
#include "physics/electromagnetism/formulation/T_Omega.h"


int main() {
    Mesh_Parser mp(Mesh_Format::GMSH);
    Mesh mesh = mp.load_mesh(SCRIPT_PATH "test_mesh_0.msh");

    T_Omega t_o(mesh);

    Logger::export_to_file("simuEM.log");
    return 0;
}