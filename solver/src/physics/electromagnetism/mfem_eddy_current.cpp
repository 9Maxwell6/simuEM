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
    Device device(device_config);
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 5. Define a finite element space on the mesh. Here we use the
    //    Raviart-Thomas finite elements of the specified order.


   //H1_FECollection HGradFEC(order, dim);
   //ND_FECollection HCurlFEC(order, dim);

    FiniteElementCollection *HGradFEC = new H1_FECollection(order, dim);
    FiniteElementCollection *HCurlFEC = new ND_FECollection(order, dim);


    FiniteElementSpace *HGrad_space = new FiniteElementSpace(mesh, HGradFEC);
    FiniteElementSpace *HCurl_space = new FiniteElementSpace(mesh, HCurlFEC);

    int Vsize_nd = HCurl_space->GetVSize();
    int Vsize_h1 = HGrad_space->GetVSize();
}