#include <iostream>
#include <fstream>
#include "Bem/basic/Bem.hpp"
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"
#include "Bem/Mesh/FittingTool.hpp"
#include "Bem/Simulation/ColocSim.hpp"

#include <cmath>
 
using namespace std;

using namespace Bem;

// this program simulates two freely oscillating bubbles. It expects
// the existence of a directory meshes in the build directory, where
// the ply outputs will be stored.

int main() {

    Mesh M;

    vector<Bem::real> vals;
    import_ply("../../first-tries/python-code/icosphere-ning/icos/ico-7.ply",M,vals);
    // the path has to be adapted such that it points to a .ply file describing a unit icosphere

    // here we prepare two different meshes M1 and M2, which we'll initialize with 
    // different initial potentials p1 and p2
    Mesh M1 = M;
    M1.translate(vec3(-2,0,0));
    PotVec p1(M.verts.size(),1.0);
    Mesh M2 = M;
    M2.translate(vec3(2,0,0));
    PotVec p2(M.verts.size(),-1.0);

    vector<PotVec> pots;
    pots.push_back(p1);
    pots.push_back(p2);

    vector<Mesh> bubbles;
    bubbles.push_back(M1);
    bubbles.push_back(M2);

    // these two functions join the two meshes and potentials together
    PotVec pot = expand_VD_to_joined(bubbles,pots);
    M = join_meshes(bubbles);

    size_t N(1500); // total number of steps
    Bem::real dp = 0.15; // time-integration precision
    
    // initialization of simulation object
    Bem::real epsilon = 10.0;
    Bem::real sigma = (epsilon-1.0)/2.0; // for equilibrium
    Bem::real gamma = 1.4;
    ColocSim sim(M,1.0,epsilon,sigma,gamma);
    sim.set_phi(pot);
    
    sim.set_damping_factor(0.5);
    Bem::real V_0(sim.get_volume()); 
    sim.set_V_0(volume(M1)); // such that V_0 in sim is the initial volume of one individual bubble (the same for both) and not the total volume

    
    size_t substeps = 1;
    for(size_t i(0);i<N;++i){
        cout << "\n----------\nCURRENT INDEX: " << i << "\n----------\n\n";

        cout << "sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << endl;
        
        
        for(size_t j(0);j<substeps;++j) {
            sim.evolve_system_RK4(dp);
            
        }

        sim.export_mesh("meshes/mesh-"+to_string(i)+".ply");

        if(i%6== 0) sim.remesh(0.2);
    }
    sim.export_mesh("meshes/mesh-"+to_string(N)+".ply");

    return 0;

}
