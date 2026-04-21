//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"


#include "math/fem/integration.h"
#include "math/fem/fem_space.h"
#include "math/fem/space_H1.h"
#include "math/fem/space_Hcurl.h"

#include "utils/logger.h"
#include "utils/util_la.h"

#include <Eigen/Dense>
#include <stdio.h>
//#include <config.h>
#include "world/mesh/mesh_parser.h"
#include "physics/electromagnetism/formulation/T_Omega.h"
#include "physics/electromagnetism/mfem_eddy_current.h"


using namespace simu;

int main(int argc, char** argv) {

    Logger::start_timer("Loading mesh");
    Mesh_Parser mp(Mesh_Format::GMSH);
    Mesh mesh = mp.load_mesh(SCRIPT_PATH "test_mesh_0_v2.2.msh");
    Logger::stop_timer("Loading mesh");

    la_kernel::initialize(&argc, &argv);

    Logger::start_timer("Initialize T-Omega solver");
    T_Omega T_O(mesh);
    Logger::stop_timer("Initialize T-Omega solver");

    

    Logger::start_timer("Assemble T-Omega matrix system");
    
    T_O.assemble_system();
    Logger::stop_timer("Assemble T-Omega matrix system");


    Logger::start_timer("Solve T-Omega matrix system");
    T_O.solve_system();
    Logger::stop_timer("Solve T-Omega matrix system");

    //MFEM_Eddy_Current(SCRIPT_PATH "test_mesh_0_v2.2.msh");

    la_kernel::finalize();


    Logger::export_to_file("simuEM.log");
    return 0;
}