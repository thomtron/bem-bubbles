#include <iostream>
#include <fstream>
#include <string>
#include "Bem/basic/Bem.hpp"
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Simulation/ColocSim.hpp"

#include <cmath>
#include <chrono>
 
using namespace std;
using namespace chrono;

using namespace Bem;

Bem::real K,Omega,Pa;

Bem::real waveform(vec3 x,Bem::real t) {
    return Pa*sin(K*x.x*0.0 - Omega*t);
}

int main(int argc, char *argv[]) {

    // This code simulates the time evolution of a bubble in an oscillating pressure field
    // (see waveform()). The initial radius and the acustic pressure are given by the first
    // two input arguments. The third input argument provides an existing path to a folder,
    // where the output .ply files shall be stored.

    if(argc != 4) {
        cerr << "invalid number of arguments!" << endl;
        return 1;
    }

    Bem::real radius   = atof(argv[1]);
    Bem::real pressure = atof(argv[2]);
    string folder = argv[3];
    if(folder.back() == '/') folder = folder.substr(0,folder.size()-1);
    folder += '/';
    cout << "radius:   " << radius << endl;
    cout << "pressure: " << pressure << endl;
    cout << "folder:   " << folder << endl;

    // physical parameters in SI units

    Bem::real p_infty = 101325.0; // N/m^2  ambient pressure                     Note: 101325.0 Pa = 1 atm by definition (see wikipedia)
    Bem::real sigma = 0.07275;    // N/m    surface tension                      @ 20°C https://de.wikipedia.org/wiki/Oberfl%C3%A4chenspannung
    Bem::real r0 = radius;         // m      initial radius = reference length    from Versluis_2010
    Bem::real c = 1481.0;         // m/s    sound speed of water                 @ 20°C https://en.wikipedia.org/wiki/Speed_of_sound  -  no good source...
    Bem::real f = 130000;         // Hz     acoustic frequency                   from Versluis_2010
    Bem::real l = c/f;            // m      acoustic wavelength
    Bem::real pa = pressure;      // N/m^2  acoustic pressure amplitude          from Versluis_2010
    Bem::real p_vap = 2300.0;     // N/m^2  water vapour pressure                @ 20°C https://en.wikipedia.org/wiki/Vapor_pressure
    Bem::real rho = 998.2067;     // kg/m^3 water density                        @ 20°C https://de.wikipedia.org/wiki/Eigenschaften_des_Wassers

    Bem::real p_ref = p_infty - p_vap;
    Bem::real t_ref = r0*sqrt(rho/p_ref);

    // dimensionless parameters:

    Bem::real P_ref = 1.0;                // reference pressure
    Bem::real Sigma = sigma/(r0*p_ref);   // surface tension
    Bem::real P_gas0 = P_ref + 2.0*Sigma; // chosen to be initially in equilibrium (R0 = 1.0)
              Omega = 2.0*M_PI*f*t_ref;   // pulsation
    Bem::real L = l/r0;                   // wave length
              K = 2.0*M_PI/L;             // wave number
    Bem::real Gamma = 7.0/5.0;            // air is a dominantly diatomic gas
              Pa = pa/p_ref;              // acoustic pressure amplitude

    Bem::real duration_max = 5.0;         // seconds
    
    cout << "P_ref =      " << P_ref << endl;
    cout << "Sigma =      " << Sigma << endl;
    cout << "P_gas0 =     " << P_gas0 << endl;
    cout << "Omega =      " << Omega << endl;
    cout << "Wavenumber = " << K << endl;
    cout << "Wavelength = " << L << endl;
    cout << "Pa =         " << Pa << endl;
    cout << "dt-min =     " << 0.1*M_PI/Omega << endl;

    cout << "time at first extremum: " << 0.5*M_PI/Omega << endl;

    // preparing simulation

    Mesh M;
    import_ply("../python_utils/icosphere/ico-10.ply",M);

    Bem::real dp = 0.02;
    size_t N(1000);
    
    ColocSim sim(M,P_ref,P_gas0,Sigma,Gamma,&waveform);
    sim.set_min_dt(0.1*M_PI/Omega);
    sim.set_phi(0.0);
    sim.set_damping_factor(0.2);
    sim.set_minimum_element_size(0.1);
    sim.set_maximum_element_size(0.9);
    Bem::real V_0(sim.get_volume());

    ofstream output(folder+"times.csv");

    size_t substeps = 4;
    for(size_t i(0);i<N;++i){
        cout << "\n----------\nCURRENT INDEX: " << i << "\n----------\n\n";
        
        cout << "sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << ", # elements: " << sim.mesh.verts.size()<< endl;
        
        output << i << ';' << sim.get_time() << ';' << sim.get_time()*t_ref << endl;
        sim.export_mesh(folder+"mesh-"+to_string(i)+".ply");

        if(i%1 == 0) sim.remesh(0.12);
        
        auto start = high_resolution_clock::now();

        for(size_t j(0);j<substeps;++j)
            sim.evolve_system(dp);

        auto end = high_resolution_clock::now();
        duration<Bem::real> dur(end-start);
        Bem::real duration = dur.count();
        /*
        if(duration > duration_max) {
            cout << "duration limit for one iteration surpassed! - ending simulation." << endl;
            N = i+1;
            break;
        }
        */
        
    }
    output << N << ';' << sim.get_time() << ';' << sim.get_time()*t_ref << endl;
    sim.export_mesh(folder+"mesh-"+to_string(N)+".ply");


    output.close();


    return 0;

}
