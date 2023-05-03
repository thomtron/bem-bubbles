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
    to_centerofmass(m);

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
    project_and_interpolate(my_uv,my_uv.verts,vals1,m,vals2);

    // filling up a .csv table with the distances of the points to the origin:
    ofstream output(filename_raw + ".csv");
    // compute mean radius (by integral):
    Bem::real mean_radius = 0.0;
    Bem::real integral = 0.0;
    for(size_t i(0);i<n_phi;++i) {
        for(size_t j(0);j<n_theta;++j) {
            Bem::real theta = (static_cast<Bem::real>(j)+0.5)*dt;
            mean_radius += my_uv.verts[i*n_theta + j].norm()*sin(theta)*dp*dt;
            integral += 2.0*sin(theta)*dp*dt;
        }
    }
    mean_radius = mean_radius/(4.0*M_PI);

    // write the .csv table to the output
    for(size_t i(0);i<n_phi;++i) {
        for(size_t j(0);j<n_theta-1;++j) {
            output << (my_uv.verts[i*n_theta + j].norm()-mean_radius)/mean_radius << ';';
        }
        output << (my_uv.verts[i*n_theta + n_theta -1].norm()-mean_radius)/mean_radius;
        output << endl;
    }
    output.close();
    return 0;
}