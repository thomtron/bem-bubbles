#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>

#include "Bem/basic/Bem.hpp"
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"
#include "Bem/Simulation/ColocSimPin.hpp"

using namespace std;
using namespace chrono;
using namespace Bem;

#if MIRROR_MESH

Bem::real K,Omega,Pa;

Bem::real waveform(vec3 x,Bem::real t) {
    return Pa*sin(K*x.x*0.0 - Omega*t);
}

// THINK ABOUT THE MIRRORING OF THE NODES AT X=0 !!
// VOLUME IST NOCH FALSCH !!

int main() {
    Mesh M;
    import_ply("semi-sphere.ply",M);

    // put all vertices with x-coord ~ 0 at the end of the mesh
    vector<size_t> permutation(M.verts.size());
    vector<size_t> inverse_permutation(M.verts.size());
    iota(permutation.begin(),permutation.end(),0);
    vector<vec3>& pos(M.verts);
    size_t k(pos.size());
    size_t npin(0);
    for(size_t i(0);i<k;++i) {
        if(abs(pos[i].x)<1e-5) {
            swap(permutation[i],permutation[k-1]);
            swap(pos[i],pos[k-1]);
            i--;
            k--;
            npin++;
        }
    }
    for(size_t i(0);i<permutation.size();++i) {
        inverse_permutation[permutation[i]] = i;
    }
    for(size_t i(0);i<M.verts.size();++i) {
        if(i>=M.verts.size()-npin) cout << "pinned: ";
        cout << M.verts[i].x << endl;
    }
    for(size_t i(0);i<M.trigs.size();++i) {
        Triplet trig(M.trigs[i]);
        trig.a = inverse_permutation[trig.a];
        trig.b = inverse_permutation[trig.b];
        trig.c = inverse_permutation[trig.c];
        M.trigs[i] = trig;
    }

    export_ply("test.ply",M);

    // perfect!

    // now lets initiate a ColocSimPin instance with this mesh

    // following lines borrowed from oscillations.cpp

    Bem::real radius   = 50e-6; // m
    Bem::real pressure = 8e4; // Pa
    string folder = "pinned-bubble-res/";
    
    cout << "radius:   " << radius << endl;
    cout << "pressure: " << pressure << endl;
    cout << "folder:   " << folder << endl;

    // physical parameters in SI units

    Bem::real p_infty = 101325.0; // N/m^2  ambient pressure                     Note: 101325.0 Pa = 1 atm by definition (see wikipedia)
    Bem::real sigma = 0.07275;    // N/m    surface tension                      @ 20°C https://de.wikipedia.org/wiki/Oberfl%C3%A4chenspannung
    Bem::real r0 = radius;        // m      initial radius = reference length    from Versluis_2010
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

    
    
    Bem::real dp = 0.005;
    size_t N(1000);
    
    ColocSimPin sim(M,npin,P_ref,P_gas0,Sigma,Gamma,&waveform);
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

        if(i%10 == 0){
            vector<size_t> inds(sim.mesh.verts.size() - npin);
            iota(inds.begin(),inds.end(),0);
            PotVec pot = sim.get_phi();
            sim.mesh = l2smooth(sim.mesh,pot,inds);
            sim.set_phi(pot);
            cout << "smoothed" << endl;
        }
        
        auto start = high_resolution_clock::now();

        for(size_t j(0);j<substeps;++j)
            sim.evolve_system_RK4(dp);

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

#else

int main() {
    cerr << "Attempted to simulate pinned bubble without mirror-kernel. please define MIRROR_MESH true. Done nothing." << endl;
    return 1;
}

#endif