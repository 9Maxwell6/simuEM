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
#include "physics/electromagnetism/mfem_eddy_current.h"

using namespace simu;

int main(int argc, char** argv) {

    Logger::start_timer("Loading mesh");
    Mesh_Parser mp(Mesh_Format::GMSH);
    Mesh mesh = mp.load_mesh(SCRIPT_PATH "test_mesh_0_v2.2.msh");
    Logger::stop_timer("Loading mesh");


    Logger::start_timer("Initialize T-Omega solver");
    T_Omega T_O(mesh);
    Logger::stop_timer("Initialize T-Omega solver");

    

    Logger::start_timer("Assemble T-Omega matrix system");
    #ifdef LOAD_PETSC
        PetscCall(PetscInitialize(&argc, &argv, 0, 0));
    #endif
    
    T_O.assemble_system();
    Logger::stop_timer("Assemble T-Omega matrix system");



    //MFEM_Eddy_Current(SCRIPT_PATH "test_mesh_0_v2.2.msh");


    #ifdef LOAD_PETSC
        PetscCall(PetscFinalize());
    #endif


    Logger::export_to_file("simuEM.log");
    return 0;
}