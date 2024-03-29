    cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
    project(3dgnome LANGUAGES C CXX CUDA)

  set(SOURCES
    src/main.cpp
    src/Anchor.cpp      
    src/ChromosomesSet.cpp  
    src/HierarchicalChromosome.cu  
    src/LooperSolver.cpp
    src/BedRegion.cpp   
    src/Cluster.cpp         
    src/InteractionArc.cpp         
    src/ParallelMarkArcs.cu
    src/BedRegions.cpp  
    src/Density.cpp         
    src/InteractionArcs.cpp        
    src/ParallelMonteCarloHeatmap.cu
    src/Chromosome.cpp  
    src/Heatmap.cpp         
    src/Settings.cpp
    thirdparty/common.cpp
    thirdparty/ini.c
    thirdparty/INIReader.cpp
    thirdparty/mtxlib.cpp
    thirdparty/rmsd.cpp
  )

  set(HEADERS
    include/Anchor.h      
    include/Chromosome.h      
    include/Density.h                 
    include/InteractionArc.h   
    include/Settings.h
    include/BedRegion.h   
    include/ChromosomesSet.h  
    include/Heatmap.h                 
    include/InteractionArcs.h
    include/BedRegions.h  
    include/Cluster.h         
    include/HierarchicalChromosome.h  
    include/LooperSolver.h
    thirdparty/common.h
    thirdparty/ini.h
    thirdparty/INIReader.h
    thirdparty/mtxlib.h
    thirdparty/rmsd.h
  )

    add_executable(cudaMMC)
    target_compile_options(cudaMMC PUBLIC -Wall)
    target_compile_features(cudaMMC PUBLIC cxx_std_17)
    set_target_properties(cudaMMC PROPERTIES CXX_STANDARD 17)
    set_target_properties(cudaMMC PROPERTIES CMAKE_CUDA_STANDARD 11)
    set_target_properties(cudaMMC PROPERTIES CMAKE_CUDA_STANDARD_REQUIRED ON)
    
    if(NOT CUDA_ARCH)
      message( SEND_ERROR "YOU NEED TO SPECIFY CUDA ARCHITECTURE E.G. -DCUDA_ARCH=\"80\" " )
    endif()

    set_target_properties(cudaMMC PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH}")

    set_target_properties(cudaMMC
                      PROPERTIES
                        CUDA_RUNTIME_LIBRARY static
                      )

    target_include_directories(cudaMMC PUBLIC 
      ${CMAKE_SOURCE_DIR}/thirdparty/
      ${CMAKE_SOURCE_DIR}/include/
    )

    target_sources(cudaMMC PUBLIC ${SOURCES} ${HEADERS})
     
    set_target_properties( cudaMMC
                           PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
     
    set_property(TARGET cudaMMC 
                 PROPERTY CUDA_SEPARABLE_COMPILATION ON)

               #target_link_libraries(particle_test PRIVATE particles)
