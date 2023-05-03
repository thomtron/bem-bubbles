#include <iostream>
#include <fstream>
#include <vector>
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Simulation/ColocSim.hpp"

 
using namespace std;

using namespace Bem;


int main() {

    // This code computes the potential u exterior of the bubble. We assume that meshes/mesh-113.ply
    // is generated with main.cpp. If changes to the parameters in main.cpp are conducted, they have
    // to be implemented in this file too (see below). The result is stored in the two files 
    // pot-ext.csv and pot_t-ext.csv (the potential and its time derivative). They can be visualized
    // using python_utils/pressure/pressure-pot-plot.py.

    Mesh M;
    vector<Bem::real> phi,psi;

    import_ply("meshes/mesh-113.ply",M,phi,psi);

    Bem::real dt = -1e-5;

    ColocSim sim(M);
    sim.set_phi(phi);
    sim.evolve_system(0.0,true);  // just recompute psi

    // simulation parameters - bust be the same as in simulation script!
    Bem::real epsilon = 10.0;
    Bem::real sigma = (epsilon-1.0)/2.0;
    Bem::real gamma = 1.4;
    Mesh M0;
    import_ply("../python_utils/icosphere/ico-7.ply",M0); // for V_0 (for time evolution)

    ColocSim simA(M,1.0,epsilon,sigma,gamma);
    simA.set_V_0(volume(M0));
    simA.set_phi(phi);

    // alternatively, one can take two simulation outputs and use their time
    // difference as dt -> one does not need to recompute a time step, but 
    // often we used a rather large time step thus we decided to recompute 
    // a timestep here.

    simA.evolve_system(dt,true);  // evolve in time for finite differences
    simA.evolve_system(0.0,true); // just recompute psi

    vector<vec3> x;
    
    for(double a(-4.0);a<=4.0001;a+=0.01){     
        for(double b(-1.5);b<=1.5001;b+=0.01){
            x.push_back(vec3(a,0.0,b));
        }
    }
    

    vector<Bem::real> phi_ext,phi_ext_A,phi_t_ext;

    // compute the potential at the positions given by x
    phi_ext_A = compute_exterior_pot(x,simA.mesh,simA.get_phi(),simA.get_psi());
    phi_ext   = compute_exterior_pot(x,sim.mesh,sim.get_phi(),sim.get_psi());
    
    // compute the time derivative of the potential using finite difference
    for(size_t i(0);i<phi_ext_A.size();++i)
        phi_t_ext.push_back((phi_ext_A[i] - phi_ext[i])/dt);


    // writing the vectors to two separate files 
    // (one for potentia and one for its time derivative)

    ofstream output("pot-ext.csv");
    for(size_t i(0);i<phi_ext.size();++i) {
        output << x[i].x << ';' << x[i].y << ';' << x[i].z << ';' << phi_ext[i] << endl;
    }
    output.close();

    output.open("pot_t-ext.csv");
    for(size_t i(0);i<phi_t_ext.size();++i) {
        output << x[i].x << ';' << x[i].y << ';' << x[i].z << ';' << phi_t_ext[i] << endl;
    }
    output.close();

    return 0;
}
