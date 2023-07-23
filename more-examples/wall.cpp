#include <iostream>
#include <fstream>
#include "../Bem/Mesh/Mesh.hpp"
#include "../Bem/Mesh/MeshIO.hpp"
#include "../Bem/Mesh/MeshManip.hpp"
#include "../Bem/Mesh/FittingTool.hpp"
#include "../Bem/Simulation/ColocSim.hpp"
#include "../Bem/Simulation/GalerkinSim.hpp"

#include <cmath>
#include <numeric>
 
using namespace std;

using namespace Bem;

int main() {

    Mesh M;

    // importing/adding meshes to M

    vector<Bem::real> vals;
    import_ply("../python_utils/icosphere/ico-10.ply",M,vals);

    Bem::real RM = 1.0; // r-max
    Bem::real R0 = 1.5142e-2;
    Bem::real gamma = 1.0;

    Bem::real DP = 1.0; // = p_inf - p_vap here
    /*
    M.scale(1.0/1.0);
    ColocSim s(M,0.0,0.0,0.0,0.0,0.0);
    s.set_phi(-1.0);
    s.evolve_system(0.0);
    return 0;
    */
    
    M.scale(R0);
    M.translate(vec3(-1.0,0.0,0.0)*gamma);

    Bem::real V_0 = volume(M);


    Bem::real dp = 0.02; //0.005;
    size_t N(90);

    
    // parameters of Wang_2014
    Bem::real epsilon = 0.0; //20.0;
    Bem::real sigma = 0.0; //0.1639103347599146; //(epsilon-1.0)/2.0; // for equilibrium
    Bem::real lambda = 0.0; //1.667; //3.0/2.0;
    ColocSim sim(M,DP,epsilon,sigma,lambda);
    sim.set_phi(-R0*sqrt(2.0/3.0*DP*(pow(RM/R0,3)-1)));

    
    size_t substeps = 10;
    for(size_t i(0);i<N;++i){
        cout << "\n----------\nCURRENT INDEX: " << i << "\n----------\n\n";
        
        cout << "sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << endl;

        
        for(size_t j(0);j<substeps;++j)
            sim.evolve_system_RK4(dp);

        sim.export_mesh("meshes-wall/mesh-"+to_string(i)+".ply");

        //if(i%10 == 0) sim.remesh(0.25);
    }
    sim.export_mesh("meshes/mesh-"+to_string(N)+".ply");


    return 0;

}