#include "math/fem/fem_space.h"


using namespace simu;


FEM_Space::FEM_Space(int dim, int p): dim_(dim), p_(p)
{
    
}


/*
FEM_Space::FEM_Space(int dim, int p, int n_node, int n_edge, int n_face, int n_volume, 
    int n_dof, int n_dof_per_node, int n_dof_per_edge, int n_dof_per_face, int n_dof_per_volume): dim_(dim), p_(p),
                                                                                                  n_node_(n_node),
                                                                                                  n_edge_(n_edge),
                                                                                                  n_face_(n_face),
                                                                                                  n_volume_(n_volume),
                                                                                                  n_dof_(n_dof),
                                                                                                  n_dof_per_node_(n_dof_per_node),
                                                                                                  n_dof_per_edge_(n_dof_per_edge),
                                                                                                  n_dof_per_face_(n_dof_per_face),
                                                                                                  n_dof_per_volume_(n_dof_per_volume){

}
*/
