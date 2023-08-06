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

add_executable(pot-ext pot-ext.cpp)
target_link_libraries(pot-ext simulation integration mesh)

add_executable(color color.cpp)
target_link_libraries(color mesh)

add_executable(pinned-bubble pinned-bubble.cpp)
target_link_libraries(pinned-bubble simulation integration mesh)

add_executable(open-mesh-test open-mesh-test.cpp)
target_link_libraries(open-mesh-test integration mesh)

add_executable(wall more-examples/wall.cpp)
target_link_libraries(wall simulation integration mesh)

add_executable(wall_check more-examples/wall_check.cpp)
target_link_libraries(wall_check simulation integration mesh)

add_executable(taib-cloud more-examples/taib-cloud.cpp)
target_link_libraries(taib-cloud simulation integration mesh)

add_executable(time-interpolation time-interpolation.cpp)
target_link_libraries(time-interpolation simulation integration mesh)

add_executable(oscillations oscillations.cpp)
target_link_libraries(oscillations simulation integration mesh)

add_executable(remeshing-for-spherical-harmonics remeshing-for-spherical-harmonics.cpp)
target_link_libraries(remeshing-for-spherical-harmonics mesh)

include_directories(${EIGEN_INCLUDE})
