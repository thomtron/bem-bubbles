This repository contains a code for simulating the evolution of bubbles in a liquid 
using the Boundary Element Method (BEM). The directory Bem contains four subdirectories
named basic, Mesh, Integration and Simulation. These direcotries contain all relevant 
header and source files and some cmake files for compiling the source files as libraries.
The main CMakeLists.txt in the root directory links these libraries to the executable of 
main.cpp. The code depends on the Eigen library. It must be of version 3.4 or higher and
it can be included by simply pasting the path to it into the file eigen_include.dir. 

The main.cpp file, which is compiled when executing 'cmake ..' and then 'make' in the 
build directory, initializes a simulation with two bubbles with opposite initial 
radial velocities. You can run it by first executing the command in python_utils/icosphere:
python ply-icosphere.py 7 ico-7.ply
which produces a .ply with seven subdividions. This file will be read by main. Then you
should go into the directory build. There you create a directory 'meshes' by:
mkdir meshes
and then you can run the program with ./main. There should be written a .ply file into 
meshes for each iteration of the simulation.

The c++ file pot-ext.cpp and the python script python_utils/pressure/pressure-pot-plot.py
can be used to visualize the flow potential and the pressure field outside the bubble.

oscillations.cpp simulates the time evolution of a bubble in an oscillating pressure field
(see waveform()). The initial radius and the acustic pressure are given by the first two 
input arguments. The third input argument provides an existing path to a folder, where the 
output .ply files shall be stored.

oscillations-array.sh can be executed with bash while being in a python environment allowing 
to import numpy, pyshtools and multiprocessing. This shell script then automatically runs
oscillations.cpp for a range of initial radii and acustic pressures and computes the spherical
harmonics power spectrum associated to the bubble surfaces. The results can be viewed using
python_utils/spherical-harmonics/visualising-sh-coeffs.py.

We further provide basic means of visualization of the .ply files in python_utils/visualization.
For visualizing the .ply meshes in a folder and following the format "mesh-#.ply", change the
symbolic link "models" in the folder python_utils/visualization/interface to a link which points
to your folder (by default it points to build/meshes). Then run the python script 
python_utils/visualization/visualization.py in a python environment where the libraries os, glob
and eel are installed.
By clicking the left mouse button you can rotate the camera around the center of 
the current view and by clicking the right mouse button you can move the center of the 
current view. The top slider allows to scroll through the meshes and the two black-bar-sliders
allow to change the color mapping of the mesh's vertex colors. If many large meshes have to be 
loaded, you have to wait a little bit until the visualization works correctly. If the interval
[0,1000] is not adequate, change the values of min and max in the input element of line 13
in index.html. This visualization tool is rather basic and meant for quick testing purposes, so
do not judge it too harshly ;)