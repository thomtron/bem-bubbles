cmake_minimum_required(VERSION 3.16)

project(BEM)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic)
set(CMAKE_BUILD_TYPE Release)
#add_compile_definitions(VERBOSE)


set(DIRECTORIES Bem Bem/Mesh Bem/Simulation Bem/Integration)

file(STRINGS eigen_include.dir EIGEN_INCLUDE LIMIT_COUNT 1)
#set(EIGEN_INCLUDE /usr/local/include/eigen-3.4.0)

message("EIGEN_INCLUDE = ${EIGEN_INCLUDE}")

foreach(directory ${DIRECTORIES})
  add_subdirectory(${directory})
endforeach(directory)


add_executable(main main.cpp)
target_link_libraries(main simulation integration mesh)

include_directories(${EIGEN_INCLUDE})




# seen in the scisim course, maybe consider looking into that...

# This needs to be in both the top-level and the test CMake file.
#enable_testing()
