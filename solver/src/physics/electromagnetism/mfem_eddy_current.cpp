#include "physics/electromagnetism/mfem_eddy_current.h"

MFEM_Eddy_Current::MFEM_Eddy_Current(const char * mesh_file, 
                                     int order, 
                                     bool pa, 
                                     const char * device_config, 
                                     bool visualization): 
                                     mesh_file_(mesh_file),
                                     order_(order),pa_(pa),
                                     device_config_(device_config),
                                     visualization_(visualization)
      
{
    mfem::Device device(device_config);
    mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 5. Define a finite element space on the mesh. Here we use the
    //    Raviart-Thomas finite elements of the specified order.


   //H1_FECollection HGradFEC(order, dim);
   //ND_FECollection HCurlFEC(order, dim);

    mfem::FiniteElementCollection *HGradFEC = new mfem::H1_FECollection(order, dim);
    mfem::FiniteElementCollection *HCurlFEC = new mfem::ND_FECollection(order, dim);


    mfem::FiniteElementSpace *HGrad_space = new mfem::FiniteElementSpace(mesh, HGradFEC);
    mfem::FiniteElementSpace *HCurl_space = new mfem::FiniteElementSpace(mesh, HCurlFEC);

    int Vsize_nd = HCurl_space->GetVSize();
    int Vsize_h1 = HGrad_space->GetVSize();
}