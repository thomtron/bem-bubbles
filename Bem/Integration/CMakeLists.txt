add_library(integration STATIC Integrator.cpp ResultTypes.cpp)

target_include_directories(integration PUBLIC ${PROJECT_SOURCE_DIR}/Bem/Integration)

include_directories(${EIGEN_INCLUDE})
