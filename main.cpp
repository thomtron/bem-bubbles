#include <iostream>
#include <fstream>
#include "Bem/basic/Bem.hpp"
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"
#include "Bem/Mesh/FittingTool.hpp"
#include "Bem/Simulation/ColocSim.hpp"
#include "Bem/Simulation/GalerkinSim.hpp"
#include "Bem/Simulation/ConConGalerkinSim.hpp"
#include "Bem/Simulation/ConLinGalerkinSim.hpp"

#include <cmath>
 
using namespace std;

using namespace Bem;

int main() {

    //Icosphere ico(3);

    Mesh M;

    // importing/adding meshes to M

    vector<Bem::real> vals;
    import_ply("../../first-tries/python-code/icosphere-ning/icos/ico-7.ply",M,vals);

    M.add(M,vec3(2.0,0.0,0.0));
    M.add(M,vec3(1.0,1.0,0.0));


    //Bem::real L = 8.2230;
    Bem::real T = 30.0;
    //size_t N(1500);
    //Bem::real dt = L/static_cast<Bem::real>(N);
    Bem::real dp = 0.04; //good value
    //Bem::real dt = 0.0025;
    //size_t N = T/dt;
    
    // parameters of Wang_2014
    Bem::real epsilon = 20.0;
    Bem::real sigma = 0.1639103347599146; //(epsilon-1.0)/2.0; // for equilibrium
    Bem::real gamma = 1.667; //3.0/2.0;
    ColocSim sim(M,1.0,epsilon,sigma,gamma);
    sim.set_phi(0.0);


    ConConGalerkinSim concon(M,1.0,epsilon,sigma,gamma);
    ConLinGalerkinSim conlin(M,1.0,epsilon,sigma,gamma);
    GalerkinSim galerkin(M,1.0,epsilon,sigma,gamma);

    concon.evolve_system(0.1);
    conlin.evolve_system(0.1);
    galerkin.evolve_system(0.1);

    //sim.set_dp_balance(6.0);
    Bem::real V_0(sim.get_volume());

    //ofstream output("times.csv");
    ofstream output("bubble-info.csv");

    
    size_t substeps = 1;
    //for(size_t i(0);i<N;++i){
    size_t i(0);
    //int remesh_ind(0);
    while(sim.get_time()<T and i<10000){
        cout << "\n----------\nCURRENT INDEX: " << i << "\n----------\n\n";

        //output << i << ';' << sim.get_time() << ';' << sim.get_volume()/V_0 << endl; // *sqrt(p_inf-p_vap)/RM
        
        cout << "sim-time: " << sim.get_time() << ", volume/V_0: " << sim.get_volume()/V_0 << endl;



        /* TODO's:
            - make good plot of taib-situation !
            - better curvature vector
            - dp-dx-plot
            - evt. multistep method
        */

        
        
        for(size_t j(0);j<substeps;++j) {
            sim.evolve_system(dp);
            //sim.remesh(0.15);
            /*
            Mesh tmp = sim.mesh;
            relax_vertices(tmp);
            vector<Bem::real> vals;
            project_and_interpolate(tmp,vals,sim.mesh,sim.get_phi());
            sim.mesh = tmp;
            sim.set_phi(vals);
            */

            vector<Bem::real> phi(sim.get_phi());
            vector<Bem::real> psi(sim.get_psi());

            Bem::real mean_rad(0.0);
            Bem::real max_rad(sim.mesh.verts[0].norm());
            Bem::real min_rad = max_rad;

            Bem::real mean_phi(0.0);
            Bem::real max_phi(phi[0]);
            Bem::real min_phi(phi[0]);

            Bem::real mean_psi(0.0);
            Bem::real max_psi(psi[0]);
            Bem::real min_psi(psi[0]);
            for(size_t i(0);i<sim.mesh.verts.size();++i) {
                Bem::real rad = sim.mesh.verts[i].norm();
                mean_rad += rad;
                min_rad = min(rad,min_rad);
                max_rad = max(rad,max_rad);
                mean_psi += psi[i];
                min_psi = min(min_psi,psi[i]);
                max_psi = max(max_psi,psi[i]);
                mean_phi += phi[i];
                min_phi = min(min_phi,phi[i]);
                max_phi = max(max_phi,phi[i]);
            }
            mean_rad /= static_cast<Bem::real>(sim.mesh.verts.size());
            mean_phi /= static_cast<Bem::real>(sim.mesh.verts.size());
            mean_psi /= static_cast<Bem::real>(sim.mesh.verts.size());
            output << i << ';' << sim.get_time() << ';' << mean_rad << ';' << min_rad << ';' << max_rad;
            output << ';' << mean_phi << ';' << min_phi << ';' << max_phi;
            output << ';' << mean_psi << ';' << min_psi << ';' << max_psi << ';' << sim.mesh.verts.size() << endl;
            
        }

        sim.export_mesh("meshes/mesh-"+to_string(i)+".ply");

        if(i%10== 0) sim.remesh(0.2);

        /*
        if(int(sim.get_time()/0.05)>remesh_ind) {
            sim.remesh(0.1);
            remesh_ind = int(sim.get_time()/0.05);
            cout << "remeshing # " << remesh_ind << endl;
        }*/

        i++;
    }
    sim.export_mesh("meshes/mesh-"+to_string(i)+".ply"); // N instead of i here normally


    output.close();


    return 0;

}
