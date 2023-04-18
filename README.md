This repository contains a code for simulating the evolution of bubbles in a liquid 
using the Boundary Element Method (BEM). The directory Bem contains four subdirectories
named basic, Mesh, Integration and Simulation. These direcotries contain all relevant 
header and source files and some cmake files for compiling the source files as libraries.
The main CMakeLists.txt in the root directory links these libraries to the executable of 
main.cpp. The code depends on the Eigen library. It must be of version 3.4 or higher and
it can be included by simply pasting the path to it into the file eigen_include.dir. 

The main.cpp file, which is compiled when executing 'cmake ..' and then 'make' in the 
build directory, initializes a simulation with two bubbles with opposite initial 
radial velocities. You can run it by first executing the command:
python ply-icosphere.py 7 ico-7.ply
which produces a .ply with seven subdividions. This file will be read by main. Then you
should go into the directory build. There you create a directory 'meshes' by:
mkdir meshes
and then you can run the program with ./main. There should be written a .ply file into 
meshes for each iteration of the simulation.
