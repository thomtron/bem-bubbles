#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <chrono>
#include <algorithm> // for clamp

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

Bem::real smoothstep(Bem::real edge0,Bem::real edge1, Bem::real x) {
    x = clamp((x-edge0)/(edge1-edge0),0.0,1.0);
    return x*x*(3.0 - 2.0*x);
}

vector<Bem::real> load_scalar_vector(string const& fname) {
    vector<Bem::real> data;
    Bem::real elm;
    ifstream input(fname);
    while(input >> elm) {
        data.push_back(elm);
    }
    input.close();
    return data;
}

Bem::real interpolate(vector<Bem::real> const& t,vector<Bem::real> const& x, Bem::real const& t0) {
    if(t0 < t.front()) return x.front();
    if(t0 > t.back()) return x.back();
    size_t i(0);
    while(t[i] < t0) i++;
    Bem::real s = (t0 - t[i-1])/(t[i]-t[i-1]);
    return x[i-1] + s*(x[i]-x[i-1]);
}


int main() {
    Mesh M;
    import_ply("../init_conditions/init-pinned.ply",M);

    // loading envelope from .csv
    //vector<Bem::real> h_ampl = load_scalar_vector("../init_conditions/r_env.csv");
    //vector<Bem::real> t_ampl = load_scalar_vector("../init_conditions/t_env.csv");

    // perfect!

    // now lets initiate a ColocSimPin instance with this mesh

    // following lines borrowed from oscillations.cpp

    Bem::real radius   = 90e-6; // m 
    Bem::real pressure = 3.6e3; // Pa
    string folder = "/cluster/home/threnggli/results/f=30e3_r=90e-6_p=3.6e3_beta=0.2_rem=0.1_epsilon=1e-2_b-nonlin-0.01-smo-fine-slow-new/";
    
    cout << "radius:   " << radius << endl;
    cout << "pressure: " << pressure << endl;
    cout << "folder:   " << folder << endl;

    // physical parameters in SI units

    Bem::real p_infty = 101325.0; // N/m^2  ambient pressure                     Note: 101325.0 Pa = 1 atm by definition (see wikipedia)
    Bem::real sigma = 0.07275;    // N/m    surface tension                      @ 20°C https://de.wikipedia.org/wiki/Oberfl%C3%A4chenspannung
    Bem::real r0 = radius;        // m      initial radius = reference length    from Versluis_2010
    Bem::real c = 1481.0;         // m/s    sound speed of water                 @ 20°C https://en.wikipedia.org/wiki/Speed_of_sound  -  no good source...
    Bem::real f = 30e3;           // Hz     acoustic frequency                   from Versluis_2010
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

    Bem::real duration_max = 60.0*5;         // seconds
    
    cout << "t_ref =      " << t_ref << endl;
    cout << "P_ref =      " << P_ref << endl;
    cout << "Sigma =      " << Sigma << endl;
    cout << "P_gas0 =     " << P_gas0 << endl;
    cout << "Omega =      " << Omega << endl;
    cout << "Wavenumber = " << K << endl;
    cout << "Wavelength = " << L << endl;
    cout << "Pa =         " << Pa << endl;
    cout << "dt-min =     " << 0.1*M_PI/Omega << endl;

    cout << "time at first extremum: " << 0.5*M_PI/Omega << endl;

    
    
    Bem::real dp = 0.003;
    size_t N(10000);
    
    ColocSimPin sim(M,P_ref,P_gas0,Sigma,Gamma,&waveform);
    sim.set_min_dt(0.1*M_PI/Omega);
    sim.set_phi(0.0);
    sim.set_damping_factor(0.5);
    sim.set_minimum_element_size(0.15); // before: 0.2
    sim.set_maximum_element_size(0.9);
    Bem::real V_0(sim.get_volume());

    ofstream output(folder+"times.csv");

    Bem::real remesh_coeff = 0.1;

    sim.remesh(remesh_coeff);
    //sim.remesh(remesh_coeff); // wegen grosser auflösungsdifferenz

    Bem::real Pa_amp = Pa;

    size_t substeps = 4;
    for(size_t i(0);i<N;++i){
        //Pa = interpolate(t_ampl,h_ampl,sim.get_time()*t_ref)*Pa_amp;
        Pa = smoothstep(0.0,0.00015,sim.get_time()*t_ref)*Pa_amp;

        cout << "\n----------\nCURRENT INDEX: " << i << "\n----------\n\n";
        
        cout << "sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << ", # elements: " << sim.mesh.verts.size()<< endl;
        
        output << i << ';' << sim.get_time() << ';' << sim.get_time()*t_ref << ';' << Pa << endl;
        stringstream fname;
        fname << folder << "mesh-" << setw(6) << setfill('0') << i << ".ply";
        sim.export_mesh(fname.str());

        if(i%6 == 0){
            vector<size_t> inds(sim.mesh.verts.size() - sim.N_pin);
            iota(inds.begin(),inds.end(),0);
            PotVec pot = sim.get_phi();
            Mesh smoothed = l2smooth(sim.mesh,pot,inds);
            CoordVec c = smoothed.verts + (-1.0)*sim.mesh.verts;
            sim.mesh.verts = sim.nopenetration(1e-2,1.0,sim.mesh.verts,c);
            sim.set_phi(pot);
            cout << "smoothed" << endl;
            //if(i%12 == 0) {
                sim.remesh(remesh_coeff);
                cout << "remeshed" << endl;
            //}
        }

        //if((i-3)%12 == 0) sim.remesh(remesh_coeff); // invalidates psi- and partly phi values
        
        auto start = high_resolution_clock::now();

        for(size_t j(0);j<substeps;++j)
            sim.evolve_system_RK4(dp);

        auto end = high_resolution_clock::now();
        duration<Bem::real> dur(end-start);
        Bem::real duration = dur.count();

        if(duration < 60.0) cout << "time for this step: " << duration << " seconds" << endl;
        else                cout << "time for this step: " << duration/60.0 << " minutes" << endl;
        
        if(duration > duration_max) {
            cout << "duration limit for one iteration surpassed! - ending simulation." << endl;
            N = i+1;
            break;
        }

        if(sim.get_time()*t_ref > 0.01) {
            cout << "maximum physical simulation time reached." << endl;
            N = i+1;
            break;
        }
        
        
    }
    output << N << ';' << sim.get_time() << ';' << sim.get_time()*t_ref << ';' << Pa << endl;
    stringstream fname;
    fname << folder << "mesh-" << setw(6) << setfill('0') << N << ".ply";
    sim.export_mesh(fname.str());


    output.close();

    return 0;
}

#else

int main() {
    cerr << "Attempted to simulate pinned bubble without mirror-kernel. please define MIRROR_MESH true. Done nothing." << endl;
    return 1;
}

#endif
