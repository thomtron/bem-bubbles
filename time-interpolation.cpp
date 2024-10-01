#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> // for replace()
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Simulation/ColocSim.hpp"

 
using namespace std;

using namespace Bem;


int main(int argc, char *argv[]) {

    // Exactly two additional arguments have to be provided: The desired timestep dt (in seconds)
    // and the directory wich contains meshes named mesh-#.ply and a file containing the timepoints
    // named times.csv

    if(argc != 3) {
        cerr << "invalid number of arguments!" << endl;
        return 1;
    }

    Bem::real dt = atof(argv[1]);
    if(dt == 0.0) {
        cerr << "invalid dt value." << endl;
        return 1;
    }
    string dir = argv[2];

    if(dir.back() == '/') dir.pop_back();
    dir.push_back('/');

    vector<Bem::real> t,t_code;
    ifstream times(dir+"times.csv");
    string line;
    while(getline(times,line)) {
        replace(line.begin(),line.end(),';',' ');
        istringstream stream(line);
        Bem::real val;
        stream >> val; cout << val << ' ';
        stream >> val; cout << val << ' ';
        t_code.push_back(val);
        stream >> val; cout << val << endl;
        t.push_back(val);
    }

    cout << system(("mkdir "+dir+"interpolated").c_str());

    ofstream new_times(dir+"interpolated/times.csv");

    size_t current_index = 0;
    size_t new_index = 0;
    for(Bem::real time(0.0); time <= t.back(); time += dt) {
        while(abs(t[current_index+1]-time) <= abs(t[current_index]-time)) current_index++; 
        // <= for the case that two consecutive times are the same due to rounding (can happen
        // in extreme situations where dt nearly goes to zero, remeshing normally helps and the sim can continue)
        
        Mesh M;
        vector<Bem::real> phi,psi;

        stringstream fname1,fname2;
        fname1 << "mesh-" << setw(6) << setfill('0') << current_index << ".ply";
        fname2 << "mesh-" << setw(6) << setfill('0') << new_index << ".ply";
        cout << system(("cp "+dir+fname1.str()+' '+dir+"interpolated/"+fname2.str()).c_str());

        new_times << new_index << ';' << t_code[current_index] << ';' << t[current_index] << endl;

        //import_ply(dir+"mesh-"+to_string(current_index)+".ply",M,phi,psi);
        //ColocSim sim(M);
        //vector<vec3> grads = sim.position_t(M,phi);
        
        //vector<vec3> normals = generate_vertex_normals(M);
        //for(size_t i(0);i<normals.size();++i) normals[i]*=psi[i];

        //M.verts = M.verts + (time-t[current_index])*normals;
        //export_ply_double(dir+"interpolated/mesh-"+to_string(new_index)+".ply",M,phi,psi);

        new_index++;
    }

    cout << "done." << endl;

    new_times.close();



    return 0;
}