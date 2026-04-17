
set(Mesh mesh.cpp
         mesh_parser.cpp
         element.cpp
         e_tetrahedron.cpp
         e_triangle.cpp
         e_edge.cpp
         e__transformation.cpp)
         
target_sources(${PROJECT_NAME} PRIVATE ${Mesh})



if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE ${Mesh})
endif()


