set(Utils logger.cpp
          util_hash.cpp
          util_petsc.cpp
          util_la.cpp)
          
target_sources(simu_core PRIVATE ${Utils})

if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE ${Utils})
endif()