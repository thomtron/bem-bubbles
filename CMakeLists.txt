cmake_minimum_required(VERSION 3.16)

project(BEM)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_BUILD_TYPE Release)
#add_compile_definitions(VERBOSE) // for more output


set(DIRECTORIES Bem/Mesh Bem/Simulation Bem/Integration)

file(STRINGS eigen_include.dir EIGEN_INCLUDE LIMIT_COUNT 1)

message("EIGEN_INCLUDE = ${EIGEN_INCLUDE}")

foreach(directory ${DIRECTORIES})
  add_subdirectory(${directory})
endforeach(directory)


add_executable(main main.cpp)
target_link_libraries(main simulation integration mesh)

include_directories(${EIGEN_INCLUDE})
