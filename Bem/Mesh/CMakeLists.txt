add_library(mesh STATIC Mesh.cpp HalfedgeMesh.cpp MeshIO.cpp MeshManip.cpp FittingTool.cpp)

target_include_directories(mesh PUBLIC ${PROJECT_SOURCE_DIR}/Bem/Mesh)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mesh PUBLIC OpenMP::OpenMP_CXX)
endif()

include_directories(${EIGEN_INCLUDE})
