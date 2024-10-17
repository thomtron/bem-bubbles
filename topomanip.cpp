#include <cmath>
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"

using namespace std;
using namespace Bem;

int main() {

    string dir = "/home/thomas/Documents/ma-thesis/paper/results/f=30e3_r=120e-6_p=9.4e3_beta=0.2_rem=0.1_epsilon=1e-2_b-nonlin-0.01-smo-fine/";
    string reference = dir + "mesh-001681.ply";
    string newfile = dir + "mesh-001681_cut.ply";
    string outfile = dir + "mesh-001681_cut_fixed.ply";

    Mesh r;
    vector<real> phi_r,psi_r;
    import_ply(reference,r,phi_r,psi_r);

    Mesh n;
    vector<real> phi_n,psi_n;
    import_ply(newfile,n);

    cout << n.verts.size() << endl;

    for(vec3 const& elm_n : n.verts) {
        //cout << elm_n << endl;
        real min_dist = (r.verts[0] - elm_n).norm2();
        size_t index = 0;
        for(size_t i(1);i<r.verts.size();++i) {
            real dist = (r.verts[i] - elm_n).norm2();
            if(dist < min_dist) {
                min_dist = dist;
                index = i;
            }
        }
        if(min_dist > 1e-5) cerr << "error too large" << endl;
        phi_n.push_back(phi_r[index]);
        psi_n.push_back(psi_r[index]);
    }

    export_ply_float(outfile,n,phi_n,psi_n);
   
    return 0;
}