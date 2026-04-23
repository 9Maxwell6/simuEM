
set(Mesh mesh.cpp
         mesh_parser.cpp
         element.cpp
         e_tetrahedron.cpp
         e_triangle.cpp
         e_edge.cpp)
         
target_sources(simu_core PRIVATE ${Mesh})



