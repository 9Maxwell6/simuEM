set(Utils logger.cpp
          util_hash.cpp
          util_petsc.cpp)
          
target_sources(${PROJECT_NAME} PRIVATE ${Utils})

if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE ${Utils})
endif()