
set(Mesh mesh.cpp
         element.cpp
         e_tetrahedron.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${Mesh})



if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE ${Mesh})
endif()


