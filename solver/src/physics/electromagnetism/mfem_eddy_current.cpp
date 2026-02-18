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
#ifdef LOAD_MFEM
    mfem::Device device(device_config);
    mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();
    Logger::info("MFEM - load mesh :"+std::string(mesh_file));

    // Step 1: Find all vertices on the conductor surface
    std::set<int> surface_vertices;

    for (int f = 0; f < mesh->GetNumFaces(); f++)
    {
        int e1, e2;
        mesh->GetFaceElements(f, &e1, &e2);

        // Skip boundary faces
        if (e1 < 0 || e2 < 0) continue;

        int attr1 = mesh->GetAttribute(e1);
        int attr2 = mesh->GetAttribute(e2);

        // Interface between conductor (1) and insulator (2)
        if ((attr1 == 1 && attr2 == 2) || (attr1 == 2 && attr2 == 1))
        {
            mfem::Array<int> face_verts;
            mesh->GetFaceVertices(f, face_verts);
            for (int i = 0; i < face_verts.Size(); i++)
                surface_vertices.insert(face_verts[i]);
        }
    }

    // Step 2: Re-tag conductor elements that touch the surface as attribute 3
    for (int e = 0; e < mesh->GetNE(); e++)
    {
        if (mesh->GetAttribute(e) != 1) continue;

        mfem::Array<int> elem_verts;
        mesh->GetElementVertices(e, elem_verts);

        for (int i = 0; i < elem_verts.Size(); i++)
        {
            if (surface_vertices.count(elem_verts[i]))
            {
                mesh->SetAttribute(e, 3);
                break;
            }
        }
    }


    mfem::Array<int> insulator_attrs;
    insulator_attrs.Append(2);

    mfem::Array<int> interior_conductor_attrs;
    interior_conductor_attrs.Append(1);

    mfem::Array<int> outer_layer_attrs;
    outer_layer_attrs.Append(3);


    // --- Create submeshes ---
    mfem::SubMesh insulator_submesh = mfem::SubMesh::CreateFromDomain(*mesh, insulator_attrs);
    mfem::SubMesh interior_submesh  = mfem::SubMesh::CreateFromDomain(*mesh, interior_conductor_attrs);
    mfem::SubMesh outer_submesh     = mfem::SubMesh::CreateFromDomain(*mesh, outer_layer_attrs);



    mfem::H1_FECollection HGrad_FEC(order, dim);
    mfem::ND_FECollection HCurl_FEC(order, dim);


    // Spaces on each submesh
    // Insulator: only H1
    mfem::FiniteElementSpace HGrad_insulator(&insulator_submesh, &HGrad_FEC);

    // Interior conductor: only Hcurl
    mfem::FiniteElementSpace HCurl_interior(&interior_submesh, &HCurl_FEC);

    // Outer layer: both H1 and Hcurl
    mfem::FiniteElementSpace HGrad_outer(&outer_submesh, &HGrad_FEC);
    mfem::FiniteElementSpace HCurl_outer(&outer_submesh, &HCurl_FEC);

   

    //mfem::FiniteElementSpace *HGrad_space = new mfem::FiniteElementSpace(mesh, HGradFEC);
    //mfem::FiniteElementSpace *HCurl_space = new mfem::FiniteElementSpace(mesh, HCurlFEC);

    int Vsize_h1_i = HGrad_insulator.GetVSize();
    int Vsize_h1_co = HGrad_outer.GetVSize();
    int Vsize_nd_co = HCurl_outer.GetVSize();
    int Vsize_nd_c = HCurl_interior.GetVSize();

    Logger::info("MFEM - HGrad_space #dofs: "+std::to_string(Vsize_h1_i+Vsize_h1_co));
    Logger::info("MFEM - HCurl_space #dofs: "+std::to_string(Vsize_nd_c+Vsize_nd_co));

    //std::cout << "H1 submesh boundary attributes: ";
    //h1_submesh.bdr_attributes.Print(std::cout);

    //std::cout << "Full mesh boundary attributes: ";
    //mesh->bdr_attributes.Print(std::cout);


    //
    //   A_CC_i         |
    //           A_CC_o | A_CI_o
    // -----------------+-----------------
    //           A_IC_o | A_II_o
    //                  |         A_II_i
    //

    double mu_insulator = 1.0;  // for now
    double mu_conductor = 1.0;  // for now
    double sigma_inv = 1.0;     // for now

    mfem::ConstantCoefficient mu_insulator_coeff(mu_insulator);
    mfem::ConstantCoefficient mu_conductor_coeff(mu_conductor);
    mfem::ConstantCoefficient sigma_inv_coeff(sigma_inv);

    // A_II insulator: (mu * u, v) in H1
    mfem::BilinearForm a_insulator(&HGrad_insulator);
    a_insulator.AddDomainIntegrator(new mfem::MassIntegrator(mu_insulator_coeff));
    a_insulator.Assemble();
    a_insulator.Finalize();

    // A_CC interior: (sigma^-1 curl u, curl v) + (mu * u, v) in Hcurl
    mfem::BilinearForm a_interior(&HCurl_interior);
    a_interior.AddDomainIntegrator(new mfem::CurlCurlIntegrator(sigma_inv_coeff));
    a_interior.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mu_conductor_coeff));
    a_interior.Assemble();
    a_interior.Finalize();

    // Outer layer: H1 block (mu * u, v)
    mfem::BilinearForm a_outer_h1(&HGrad_outer);
    a_outer_h1.AddDomainIntegrator(new mfem::MassIntegrator(mu_conductor_coeff));
    a_outer_h1.Assemble();
    a_outer_h1.Finalize();

    // Outer layer: Hcurl block (sigma^-1 curl u, curl v) + (mu * u, v)
    mfem::BilinearForm a_outer_hcurl(&HCurl_outer);
    a_outer_hcurl.AddDomainIntegrator(new mfem::CurlCurlIntegrator(sigma_inv_coeff));
    a_outer_hcurl.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mu_conductor_coeff));
    a_outer_hcurl.Assemble();
    a_outer_hcurl.Finalize();

    // Outer layer: coupling (mu * grad u, v), H1 trial, Hcurl test
    mfem::MixedBilinearForm a_coupling(&HGrad_outer, &HCurl_outer);
    a_coupling.AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(mu_conductor_coeff));
    a_coupling.Assemble();
    a_coupling.Finalize();


    // Outer layer: coupling (mu * u, grad v), Hcurl trial, H1 test
    // (mu * u, grad v) = (mu * grad u, v)^T
    mfem::SparseMatrix *a_coupling_T = mfem::Transpose(a_coupling.SpMat());

    // Define block offsets
    // Block 0: HCurl_interior
    // Block 1: HCurl_outer
    // Block 2: HGrad_outer
    // Block 3: HGrad_insulator
    mfem::Array<int> block_offsets(5); // num_blocks + 1
    block_offsets[0] = 0;
    block_offsets[1] = HCurl_interior.GetTrueVSize();
    block_offsets[2] = HCurl_outer.GetTrueVSize();
    block_offsets[3] = HGrad_outer.GetTrueVSize();
    block_offsets[4] = HGrad_insulator.GetTrueVSize();
    block_offsets.PartialSum();


    //
    //   A_CC_i         |
    //           A_CC_o | A_CI_o
    // -----------------+-----------------
    //           A_IC_o | A_II_o
    //                  |         A_II_i
    //
    // A_CC_i -> a_interior
    // A_CC_o -> a_outer_hcurl
    // A_II_o -> a_outer_h1
    // A_II_i -> a_insulator
    //
    // A_CI_o -> a_coupling
    // A_IC_o -> a_coupling_T
    mfem::BlockMatrix A(block_offsets);

    // Diagonal blocks
    A.SetBlock(0, 0, &a_interior.SpMat());
    A.SetBlock(1, 1, &a_outer_hcurl.SpMat());
    A.SetBlock(2, 2, &a_outer_h1.SpMat());
    A.SetBlock(3, 3, &a_insulator.SpMat());

    // Coupling blocks (outer layer: Hcurl <-> H1)
    A.SetBlock(1, 2, &a_coupling.SpMat());
    A.SetBlock(2, 1, a_coupling_T);


    

    /*
    
    //
    //   [ A_CC   A_CI ]
    //   [ A_IC   A_II ]
    //

    mfem::BilinearForm A_CC(&HCurl_space);
    // (sigma^-1(x)*curl u, curl v)
    double sigma_value = 1.;  // for now
    mfem::ConstantCoefficient sigma_inv(1.0 / sigma_value);
    A_CC.AddDomainIntegrator(new mfem::CurlCurlIntegrator(sigma_inv));

    // (mu(x)*u, v)
    double mu_conductor = 1.0; // for now
    mfem::ConstantCoefficient mu(mu_conductor);  
    A_CC.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mu));

    A_CC.Assemble();
    A_CC.Finalize();



    // H1 space contains attribute 2 and 3,  
    //2 ->insulating region,   3 -> conductor
    double mu_insulator = 1.0; // for now
    mfem::Vector mu_vals(h1_submesh.attributes.Max());
    mu_vals = 0.0;
    mu_vals(1) = mu_insulator;   // attr 2 -> index 1 (0-based)
    mu_vals(2) = mu_conductor;   // attr 3 -> index 2 (0-based)

    mfem::PWConstCoefficient mu_pw(mu_vals);
    
    mfem::BilinearForm A_II(&HGrad_space);
    A_II.AddDomainIntegrator(new mfem::MassIntegrator(mu_pw));
    A_II.Assemble();
    A_II.Finalize();


    // for A_CI and A_IC, we need new sub_mesh:
    mfem::Array<int> overlap_attrs;
    overlap_attrs.Append(3);
    mfem::SubMesh overlap_submesh = mfem::SubMesh::CreateFromDomain(*mesh, overlap_attrs);
    mfem::FiniteElementSpace HGrad_overlap(&overlap_submesh, &HGrad_FEC);
    mfem::FiniteElementSpace HCurl_overlap(&overlap_submesh, &HCurl_FEC);
    mfem::ConstantCoefficient mu_coeff(mu_conductor);

    mfem::MixedBilinearForm a_CI_local(&HGrad_overlap, &HCurl_overlap);
    a_CI_local.AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(mu_coeff));
    a_CI_local.Assemble();
    a_CI_local.Finalize();

    */

#else
    Logger::error("MFEM_Eddy_Current - MFEM support not compiled!");
    throw std::runtime_error("MFEM support not compiled!");
#endif
        
}