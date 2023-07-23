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
    import_ply("bubblecloud-random-0_01.ply",M);

    // physical parameters in SI units

    Bem::real r0 = 100*1e-6;      // m      approx. cloud diameter
    Bem::real p_infty = 101325.0; // N/m^2  ambient pressure                     Note: 101325.0 Pa = 1 atm by definition (see wikipedia)
    Bem::real p_vap = 2300.0;     // N/m^2  water vapour pressure                @ 20°C https://en.wikipedia.org/wiki/Vapor_pressure
    Bem::real rho = 998.2067;     // kg/m^3 water density                        @ 20°C https://de.wikipedia.org/wiki/Eigenschaften_des_Wassers

    Bem::real p_ref = p_infty - p_vap;
    Bem::real t_ref = r0*sqrt(rho/p_ref);

    // parameters in code units

    Bem::real RM = 0.4; // r-max
    Bem::real R0 = 1.0e-2;
    Bem::real gamma = 1.25;

    Bem::real DP = 1.0; // = p_inf - p_vap here

    M.translate(vec3(-1.0,0.0,0.0)*gamma);


    Bem::real dp = 0.005; //0.005;
    size_t N(3000);

    
    // parameters of Wang_2014
    Bem::real epsilon = 0.0; //20.0;
    Bem::real sigma = 0.0; //0.1639103347599146; //(epsilon-1.0)/2.0; // for equilibrium
    Bem::real lambda = 0.0; //1.667; //3.0/2.0;
    ColocSim sim(M,DP,epsilon,sigma,lambda);
    sim.set_damping_factor(0.6);
    sim.set_phi(-R0*sqrt(2.0/3.0*DP*(pow(RM/R0,3)-1)));
    //sim.set_phi(vals);
    Bem::real V_0(sim.get_volume());
    

    sim.remesh(0.25);

    string folder = "taib-cloud-res/";
    cout << "creating directory: " << folder << system(("mkdir "+folder).c_str());

    ofstream output(folder+"times.csv");
    
    size_t substeps = 1;
    for(size_t i(0);i<N;++i){

        cout << "Index: " << i << ", sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << ", # elements: " << sim.mesh.verts.size()<< endl;
        
        output << i << ';' << sim.get_time() << ';' << sim.get_time()*t_ref << endl;
        
        for(size_t j(0);j<substeps;++j) {
            //sim.evolve_system_RK4(dp); // fixed dt = 0.01, can change this in ColocSim.cpp
            sim.evolve_system(dp);
        }

        sim.export_mesh(folder + "mesh-"+to_string(i)+".ply");

        if(i%20 == 0) sim.remesh(0.25);
    }

    output.close();

    return 0;

}