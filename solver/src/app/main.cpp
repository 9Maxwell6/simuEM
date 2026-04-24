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
#include <functional>

//#include <config.h>
#include "world/mesh/mesh_parser.h"
#include "physics/electromagnetism/formulation/T_Omega.h"
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
    bool l2_test_flag = false;


    const std::vector<Option> options = {
        {"mesh",    'm',  true,  [&](std::string_view v) { mesh_file = std::string(v); }},
        {"l2-test", '\0',  false, [&](std::string_view)   { l2_test_flag = true;       }},
        // add more here — one line each
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


    if(l2_test_flag)
    {
        const std::vector<std::pair<std::string, double>> mesh_sweep = {
            {"test_cc_0.geo", 0.500000000},
            {"test_cc_1.geo", 0.353553391},
            {"test_cc_2.geo", 0.250000000},
            {"test_cc_3.geo", 0.176776695},
            {"test_cc_4.geo", 0.125000000},
            {"test_cc_5.geo", 0.088388348},
        };

        const std::string dat_path = DATA_OUTPUT_DIR + std::string("/T_Omega_l2.dat");
        std::ofstream l2_convergence(dat_path);
        l2_convergence << "# h                        L2_error\n";
        l2_convergence << std::scientific << std::setprecision(15);
        

        for (const auto& [mesh_file, h] : mesh_sweep) {
            scalar_t l2_error;
            {
                Logger::start_timer("Loading mesh");
                Mesh_Parser mp(Mesh_Format::GMSH);
                Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
                Logger::stop_timer("Loading mesh");

                Logger::start_timer("Initialize CurlCurl solver");
                T_Omega T_O(mesh);
                Logger::stop_timer("Initialize CurlCurl solver");

                Logger::start_timer("Assemble CurlCurl matrix system");
                T_O.assemble_system();
                Logger::stop_timer("Assemble CurlCurl matrix system");

                Logger::start_timer("Solve CurlCurl matrix system");
                T_O.solve_system();
                Logger::stop_timer("Solve CurlCurl matrix system");

                Logger::start_timer("Compute L2 error.");
                l2_error = T_O.compute_L2_error();
                Logger::stop_timer("Compute L2 error.");
            }

            std::ostringstream ss;
            ss << std::scientific << std::setprecision(15) << l2_error;
            Logger::info("[CurlCurl] h = " + std::to_string(h) + "  L2 error: " + ss.str());

            l2_convergence << h << "  " << l2_error << "\n";
            l2_convergence.flush();   // persist after every run — a crash on the
                                // finest mesh won't lose the earlier points
        }

        l2_convergence.close();

    }else{
        Logger::start_timer("Load mesh");
        Mesh_Parser mp(Mesh_Format::GMSH);
        Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
        Logger::stop_timer("Load mesh");


        Logger::start_timer("Initialize T-Omega solver");
        T_Omega T_O(mesh);
        Logger::stop_timer("Initialize T-Omega solver");

        Logger::start_timer("Assemble T-Omega matrix system");
        T_O.assemble_system();
        Logger::stop_timer("Assemble T-Omega matrix system");

        Logger::start_timer("Solve T-Omega matrix system");
        T_O.solve_system();
        Logger::stop_timer("Solve T-Omega matrix system");

        Logger::start_timer("Compute L2 error.");
        scalar_t l2_error = T_O.compute_L2_error();
        Logger::stop_timer("Compute L2 error.");

        std::ostringstream ss;
        ss << std::scientific << std::setprecision(15) << l2_error;
        Logger::info("[T-Omega] - test case L2 error: " + ss.str());

    }

    la_kernel::finalize();
    
    Logger::export_to_file("simuEM.log");
    return 0;
}