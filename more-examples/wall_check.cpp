#include <vector>
#include "../Bem/Mesh/Mesh.hpp"
#include "../Bem/Mesh/MeshIO.hpp"
#include "../Bem/Simulation/ColocSim.hpp"

 
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

    import_ply("meshes-wall/mesh-89.ply",M,phi,psi);

    ColocSim sim(M);
    sim.set_phi(phi);
    sim.evolve_system(0.0,true); // recompute psi

    vector<vec3> points;
    Bem::real incr = 0.05;
    for(Bem::real x(0.0);x<1000.0;x += incr) {
        incr*=1.8;
        points.push_back(vec3(x,0.0,0.0));
    }

    vector<Bem::real> phi_ext = compute_exterior_pot(points,sim.mesh,sim.get_phi(),sim.get_psi());

    for(size_t i(0);i<points.size();++i) {
        cout << points[i] << " : phi = " << phi_ext[i] << endl;
    }

    return 0; 
}