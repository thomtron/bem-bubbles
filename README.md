This repository contains a code for simulating the evolution of bubbles in a liquid 
using the Boundary Element Method (BEM). The directory Bem contains four subdirectories
named basic, Mesh, Integration and Simulation. These direcotries contain all relevant 
header and source files and some cmake files for compiling the source files as libraries.
The main CMakeLists.txt in the root directory links these libraries to the executable of 
main.cpp. The code depends on the Eigen library. It must be of version 3.4 or higher and
it can be included by simply pasting the path to it into the file eigen_include.dir. 
