# Add sub-projects
add_subdirectory(material)
add_subdirectory(mesh)
add_subdirectory(structure)


target_sources(${PROJECT_NAME} PRIVATE 
    element.cpp
    tetrahedron.cpp)