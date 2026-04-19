#pragma once

#include "math/fem/block.h"
#include "math/data_format.h"

#include "utils/logger.h"
#include "utils/util_la.h"

// for test only
//#include "utils/util_petsc.h"



namespace simu {

/**
 * 
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 *      |     |     |     |
 *      +-----+-----+-----+
 */
class Block_Rack
{
    friend class  FEM_System;
    friend struct Dirichlet_BC;
private:
    
    size_t n_row_; 
    size_t n_col_; 

    std::vector<Block*>  rack_;

    std::vector<size_t> unit_row_length_;
    std::vector<size_t> unit_col_length_;

    G_Matrix lhs_;
    G_Vector rhs_;
    G_Vector   x_;    // solution

    
public:
    Block_Rack(){};
    Block_Rack(size_t n_row, size_t n_col);

    void set_grid(size_t n_row, size_t n_col);

    bool insert_block(Block& block, size_t row, size_t col);
    bool compute_block_offset();

    bool is_block_ready(const Block* block) const;

    void build_linear_system();

    void solve(){
        // for text
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);

        // Set the operator (matrix A)
        KSPSetOperators(ksp, lhs_, lhs_);

        // Set solver type to Conjugate Gradient
        KSPSetType(ksp, KSPCG);

        // Optionally configure the preconditioner (e.g., Jacobi)
        PC pc;
        KSPGetPC(ksp, &pc);    // <-- add this line
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");

        // Optionally set tolerances
        KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        // Allow command-line overrides (e.g., -ksp_type, -pc_type)
        KSPSetFromOptions(ksp);

        // Solve
        KSPSolve(ksp, rhs_, x_);

        // Check convergence
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "KSP did not converge: reason %d\n", reason);
        } else {
            PetscInt its;
            KSPGetIterationNumber(ksp, &its);
            PetscPrintf(PETSC_COMM_WORLD, "Converged in %d iterations\n", its);
        }
        petsc_util::petsc_save_ascii_mat(lhs_, "lhs_mat.txt");
        petsc_util::petsc_save_ascii_vec(x_, "x_vec.txt");
        petsc_util::petsc_save_ascii_vec(rhs_, "rhs_vec.txt");
        // Clean up
        KSPDestroy(&ksp);
    }

    void delete_data();
    

    std::string print_block_rack() const;

    //TODO:
    //void initial_guess();

    //TODO: Dirichlet BC



    



};


}