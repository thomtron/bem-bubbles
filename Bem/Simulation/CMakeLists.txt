add_library(simulation STATIC Simulation.cpp ColocSim.cpp GalerkinSim.cpp LinLinSim.cpp ConConGalerkinSim.cpp ConLinGalerkinSim.cpp)

target_include_directories(simulation PUBLIC ${PROJECT_SOURCE_DIR}/Bem/Simulation)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(simulation PUBLIC OpenMP::OpenMP_CXX)
endif()

include_directories(${EIGEN_INCLUDE})
