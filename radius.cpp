#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm> // for replace

#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"

using namespace std;
using namespace Bem;


int main(int argc, char *argv[]) {
    /*
    This script imports a .ply with phi- and psi data stored as floating point values and 
    exports a .ply with color data corresponding to the imported color map and the mapping details
    implemented below
    */

    if(argc != 2) {
        cerr << "invalid number of arguments!" << endl;
        return 1;
    }

    string folder = argv[1];
    if(folder.back() == '/') folder = folder.substr(0,folder.size()-1);
    folder += '/';

    ifstream input(folder+"times.csv");
    string line;
    while(getline(input,line)) {
        replace(line.begin(),line.end(),';',' ');
        stringstream linestream(line);
        unsigned int index(0);
        double time_code(0.0);
        double time_real(0.0);
        linestream >> index;
        linestream >> time_code;
        linestream >> time_real;

        Mesh m;
        vector<Bem::real> phi,psi;
        import_ply(folder+"mesh-"+to_string(index)+".ply",m,phi,psi);
        double radius = pow(volume(m)/(4.0/3.0*M_PI),1.0/3.0);
        cout << time_real << ';' << radius << ';' << endl;
    }

}