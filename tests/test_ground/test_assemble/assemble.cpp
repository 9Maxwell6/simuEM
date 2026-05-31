#include "assemble.h"

//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"


#include "math/fem/fem_space.h"
#include "math/fem/space_Hcurl.h"
#include "math/fem/space_Hcurl.h"

#include "utils/logger.h"
#include "utils/util_la.h"

#include <Eigen/Dense>
#include <stdio.h>
//#include <config.h>
#include "world/mesh/mesh_parser.h"

// for test
#include <fstream>

using namespace simu;

Assemble::Assemble(Mesh& mesh) : mesh_(mesh), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    Hcurl_   = Hcurl_Space(mesh.get_mesh_dimension(), 1);

    H1_v_    = H1_Space(mesh.get_mesh_dimension(), 1, true, 0);
    H1_s_    = H1_Space(mesh.get_mesh_dimension(), 1);


    key_true_boundary_ = mesh_.get_keys_true_boundary()[0];  // there must be only one true boundary.
    std::string true_boundary_description = mesh.get_group_description(key_true_boundary_);
    Logger::info("Assemble - Found simulation boundary: " + true_boundary_description);
    
    
    
    for(const Key& mesh_key : mesh_.get_keys_domain())
    {
        std::string description = mesh.get_group_description(mesh_key);

        key_interior_ = mesh_key;
        mesh_.set_group_property_id(mesh_key, Domain::INTERIOR);
        break;
    }

    
    Logger::info("[Assemble] - Assign space Hcurl to the interior.");
    dof_field_ = fe_system_.register_FE_space(H1_s_, key_interior_);
    
    Logger::info("[Assemble] - register Dirichlet BC to the true boundary.");
    bc_ = fe_system_.register_Dirichlet_BC(dof_field_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    
    Logger::block_info(dof_field_.id, dof_field_.row_offset, dof_field_.col_offset, dof_field_.row_size, dof_field_.col_size);


    fe_system_.delete_block_hash();

    br_system_ = fe_system_.initialize_block_rack(1, 1);
    
    br_system_.insert_block(dof_field_,         0, 0);

    br_system_.compute_block_offset();

    
}





bool Assemble::assemble_system()
{

    // auto& e_data - per-element data shared across integrators.
    // 
    // see assemble_data.h -> Element_Data<phy_dim, ref_dim>
    //
    // Available members:
    //   e             - const Element*, current element
    //   J             - Matrix<phy_dim, ref_dim>, Jacobian at current quad point
    //   inv_J         - Matrix<ref_dim, phy_dim>, inverse/pseudo-inverse of Jacobian
    //   det_J         - double, determinant of Jacobian
    //   b_shape       - Basis_Shape, element geometry
    //   shape_space_1 - const FEM_Space*, trial function space
    //   shape_space_2 - const FEM_Space*, test function space
    //
    //
    // Template parameters:
    //   phy_dim - physical space dimension (1, 2, or 3)
    //   ref_dim - reference element dimension (1, 2, or 3), ref_dim <= phy_dim
    //
    // Note: accessed via auto& in user callbacks. Use e_data.e, e_data.J, etc.

    
    Logger::info("[Assemble] - assemble matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_field_), [&](auto& e_data, auto& mat) {

        Integrator__s_grad_S__grad_S::assemble_element_matrix(1., e_data, mat);
    
    });



    Logger::info("[Assemble] - assemble vector.");
    assemble_vec(fe_system_.assemble_vec_data(dof_field_), [&](auto& e_data, auto& vec) {

        
    
    });

    la_kernel::destroy_mat(dof_field_.mat);
    la_kernel::destroy_vec(dof_field_.vec);



    
    
    return true;
}



using namespace simu;

int main(int argc, char** argv) {

    std::vector<char*> petsc_argv_list;
    petsc_argv_list.push_back(argv[0]);

    std::string mesh_file = "test_cube_3.msh";
    bool enable_preconditioner = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if(arg.rfind("--mesh=", 0) == 0) {
            mesh_file = arg.substr(7);
        }else{
            petsc_argv_list.push_back(argv[i]);  // leave it for PETSc
        }
    }

    char** petsc_argv = petsc_argv_list.data();
    int petsc_argc = petsc_argv_list.size();
    la_kernel::initialize(&petsc_argc, &petsc_argv);

   



    Logger::start_timer("Loading mesh");
    Mesh_Parser mp(Mesh_Format::GMSH);
    Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
    Logger::stop_timer("Loading mesh");

    Logger::start_timer("Initialize assemble system");
    Assemble Assemble(mesh);
    Logger::stop_timer("Initialize assemble system");

    Logger::start_timer("Assemble system");
    Assemble.assemble_system();

    double assemble_time = Logger::stop_timer("Assemble system");

    Logger::info("[assemble time] : "+std::to_string(assemble_time));

    la_kernel::finalize();

    Logger::export_to_file("simuEM.log");
    return 0;
}
