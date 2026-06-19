//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"


#include "math/fem/fem_space.h"
#include "math/fem/space_H1.h"
#include "math/fem/space_Hcurl.h"
#include "math/data_format.h"

#include "utils/logger.h"
#include "utils/util_la.h"

#include <Eigen/Dense>
#include <stdio.h>
#include <functional>

//#include <config.h>
#include "world/mesh/mesh_parser.h"
#include "physics/electromagnetism/eddy_current/T_Omega.h"
#include "physics/electromagnetism/mfem_eddy_current.h"


struct Option 
{
    std::string_view long_name; 
    char             short_name; 
    bool             takes_value;
    std::function<void(std::string_view)> handler;
};

using namespace simu;

int main(int argc, char** argv) 
{
    std::string mesh_file    = "test_cc.msh";
    bool h_refine = false;
    bool l2_test = false;
    bool enable_preconditioner = false;


    const std::vector<Option> options = {
        {"mesh",     'm',   true,  [&](std::string_view v) { mesh_file = std::string(v);    }},
        {"h-refine", 'h',   false, [&](std::string_view)   { h_refine = true;               }},
        {"l2-test",  '\0',  false, [&](std::string_view)   { l2_test = true;                }},
        {"pc",       '\0',  false, [&](std::string_view)   { enable_preconditioner = true;  }},
    };

    std::vector<char*> petsc_argv_list{ argv[0] };

    for (int i = 1; i < argc; ++i) 
    {
        std::string_view arg = argv[i];
        bool matched = false;

        for (const auto& opt : options) 
        {
            std::string long_flag  = "--" + std::string(opt.long_name);
            std::string short_flag = opt.short_name ? std::string("-") + opt.short_name: std::string{};

            if (opt.takes_value) {
                if (arg.rfind(long_flag + "=", 0) == 0) {
                    opt.handler(arg.substr(long_flag.size() + 1));
                    matched = true;
                    break;
                }
                if (opt.short_name && arg.rfind(short_flag + "=", 0) == 0) {
                    opt.handler(arg.substr(short_flag.size() + 1));
                    matched = true;
                    break;
                }
            } else {
                if (arg == long_flag || (opt.short_name && arg == short_flag)) {
                    opt.handler("");
                    matched = true;
                    break;
        }}}
        if (!matched) petsc_argv_list.push_back(argv[i]);
    }

    int    petsc_argc = static_cast<int>(petsc_argv_list.size());
    char** petsc_argv = petsc_argv_list.data();
    la_kernel::initialize(&petsc_argc, &petsc_argv);



    if(h_refine)
    {
        const std::vector<std::pair<std::string, double>> mesh_sweep = {
            //{"test_cc_0.geo", 0.500000000},

            {"test_cc_1.geo", 0.353553391},
            {"test_cc_2.geo", 0.250000000},
            {"test_cc_3.geo", 0.176776695},
            {"test_cc_4.geo", 0.125000000},
            {"test_cc_5.geo", 0.088388348},
            {"test_cc_6.geo", 0.062500000},
            {"test_cc_7.geo", 0.044200000},
            {"test_cc_8.geo", 0.031250000},

            //{"test_torus_1.geo", 0.353553391},
            //{"test_torus_2.geo", 0.250000000},
            //{"test_torus_3.geo", 0.176776695},
            //{"test_torus_4.geo", 0.125000000},
            //{"test_torus_5.geo", 0.088388348},
            //{"test_torus_6.geo", 0.062500000},
            //{"test_torus_7.geo", 0.044200000},
            //{"test_torus_8.geo", 0.031250000},

            //{"test_cc_9.geo", 0.022100000},
            //{"test_cc_10.geo", 0.01562500}
        };

        std::string file_str = enable_preconditioner ? "/T_Omega_pc.dat" : "/T_Omega.dat";

        const std::string dat_path = DATA_OUTPUT_DIR + file_str;
        std::ofstream l2_convergence(dat_path);
        l2_convergence << "# h                    #cell      L2_error           #iteration      condition_number     assemble[s]            solve[s]\n";
        l2_convergence << std::scientific << std::setprecision(15);
        

        for (const auto& [mesh_file, h] : mesh_sweep) {
            petsc_util::petsc_print_memory_usage("iter N start");
            scalar_t l2_error;
            size_t n_element = 0;
            int n_iteration = -1;
            double condition_numbder = 0.;
            double assemble_time = 0.;
            double solving_time = 0.;
            {
                Logger::start_timer("Loading mesh");
                Mesh_Parser mp(Mesh_Format::GMSH);
                Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
                n_element = mesh.get_mesh_elements().size();
                Logger::stop_timer("Loading mesh");

                Logger::start_timer("Initialize T-Omega solver");
                T_Omega T_O(mesh, enable_preconditioner);
                Logger::stop_timer("Initialize T-Omega solver");

                Logger::start_timer("Assemble T-Omega matrix system");
                T_O.assemble_system();
                assemble_time += Logger::stop_timer("Assemble T-Omega matrix system");

                if(enable_preconditioner){
                    Logger::start_timer("Assemble T-Omega preconditioner");
                    T_O.assemble_preconditioner();
                    assemble_time += Logger::stop_timer("Assemble T-Omega preconditioner");
                }

                for(int i=0; i<1; ++i)
                {

                    Logger::start_timer("Solve T-Omega matrix system");
                    T_O.solve_system();
                    solving_time = Logger::stop_timer("Solve T-Omega matrix system");

                    n_iteration = T_O.get_n_iteration();
                    condition_numbder = T_O.get_system_condition();

                    if(l2_test){
                        Logger::start_timer("Compute L2 error.");
                        l2_error = T_O.compute_L2_error();
                        Logger::stop_timer("Compute L2 error.");
                    }

                    std::ostringstream ss;
                    ss << std::scientific << std::setprecision(15) << l2_error;
                    Logger::info("[T-Omega] h = " + std::to_string(h) + "  L2 error: " + ss.str());

                    l2_convergence << h << "  " << n_element 
                                        << "  " << l2_error 
                                        << "  " << n_iteration
                                        << "  " << condition_numbder
                                        << "  " << assemble_time
                                        << "  " << solving_time << "\n";
                    l2_convergence.flush();  

                }

                
                

                T_O.finalize();
            }

            

            petsc_util::petsc_print_memory_usage("iter N end");
            //PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
        }

        l2_convergence.close();

    }else{
        PetscLogDefaultBegin();
        Logger::start_timer("Load mesh");
        Mesh_Parser mp(Mesh_Format::GMSH);
        Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
        Logger::stop_timer("Load mesh");


        Logger::start_timer("Initialize T-Omega solver");
        T_Omega T_O(mesh, enable_preconditioner);
        Logger::stop_timer("Initialize T-Omega solver");

        Logger::start_timer("Assemble T-Omega matrix system");
        T_O.assemble_system();
        Logger::stop_timer("Assemble T-Omega matrix system");

        if(enable_preconditioner){
            Logger::start_timer("Assemble T-Omega preconditioner");
            T_O.assemble_preconditioner();
            Logger::stop_timer("Assemble T-Omega preconditioner");
        }

        Logger::start_timer("Solve T-Omega matrix system");
        T_O.solve_system();
        Logger::stop_timer("Solve T-Omega matrix system");

        if(l2_test){
            Logger::start_timer("Compute L2 error.");
            scalar_t l2_error = T_O.compute_L2_error();
            Logger::stop_timer("Compute L2 error.");

            std::ostringstream ss;
            ss << std::scientific << std::setprecision(15) << l2_error;
            Logger::info("[T-Omega] - test case L2 error: " + ss.str());
        }

        

        T_O.finalize();
    }

    //PetscLogView(PETSC_VIEWER_STDOUT_WORLD);

    la_kernel::finalize();
    
    Logger::export_to_file("simuEM.log");
    return 0;
}