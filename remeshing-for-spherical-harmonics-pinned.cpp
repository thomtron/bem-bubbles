#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"

using namespace std;
using namespace Bem;

int main(int argc, char *argv[]) {

    // This code creates a uv sphere, reads the .ply file
    // given in the input argument and puts its center of 
    // mass to the origin. It then projects the vertices
    // of the uv sphere on this mesh and stores the 
    // relative deviations of the distances of the vertices 
    // to the origin from their mean distance into a .csv
    // file with the same name as the .ply file.

    // version with pinned bubble

    // used in oscillations-array.sh

    if(argc != 2) {
        cerr << "invalid number of arguments!" << endl;
        return 1;
    }

    string filename = argv[1];
    size_t pos(filename.find(".ply"));
    string filename_raw = filename.substr(0,pos);
    cout << filename << ' ' << filename_raw << endl;

    Mesh m;
    import_ply(filename,m);
    m.rotate(vec3(0.0,-M_PI_2,0.0));
    vec3 com = centerofmass(m); // will be used further down...
    for(vec3& pos : m.verts)
        pos -= com;

    // create a uv-sphere with 128 phi-subdivisions and 64 theta-subdivisions
    Mesh my_uv;
    size_t n_phi(128);
    size_t n_theta(64);
    Bem::real dt = M_PI/static_cast<Bem::real>(n_theta);
    Bem::real dp = 2.0*M_PI/static_cast<Bem::real>(n_phi);
    for(size_t i(0);i<n_phi;++i) {
        for(size_t j(0);j<n_theta;++j) {
            Bem::real phi = (static_cast<Bem::real>(i)+0.5)*dp;
            Bem::real theta = (static_cast<Bem::real>(j)+0.5)*dt;
            vec3 pos(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
            my_uv.verts.push_back(pos);
        }
    }
    
    // now project this uv-sphere on the imported mesh
    vector<Bem::real> vals1(my_uv.verts.size()),vals2(m.verts.size()); // dummy variables
    // my_uv.verts are the normals to my_uv !
    project_from_origin(my_uv.verts,m,com.z);

    export_ply("prepared-mesh-sh.ply",my_uv);
    export_ply("base-mesh-sh.ply",m);

    // filling up a .csv table with the distances of the points to the origin:
    ofstream output(filename_raw + ".csv");

    // write the .csv table to the output
    for(size_t i(0);i<n_phi;++i) {
        for(size_t j(0);j<n_theta-1;++j) {
            output << my_uv.verts[i*n_theta + j].norm() << ';';
        }
        output << my_uv.verts[i*n_theta + n_theta -1].norm();
        output << endl;
    }
    output.close();
    return 0;
}