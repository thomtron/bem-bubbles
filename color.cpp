#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"

using namespace std;
using namespace Bem;


vec3 interpolate(vector<Bem::real> const& pos,vector<vec3> const& cmap,Bem::real x) {
    size_t n(pos.size());
    if(x<=pos[0]) return cmap[0];
    if(x>=pos[n-1]) return cmap[n-1];
    size_t i(1);
    while(x > pos[i]) ++i;
    Bem::real z = (x-pos[i-1])/(pos[i]-pos[i-1]);
    return cmap[i-1]*z + cmap[i]*(1-z);
}

int main(int argc, char *argv[]) {
    /*
    This script imports a .ply with phi- and psi data stored as floating point values and 
    exports a .ply with color data corresponding to the imported color map and the mapping details
    implemented below
    */

    if(argc != 4 && argc != 6) {
        cerr << "invalid number of arguments!" << endl;
        return 1;
    }

    string model_fname = argv[1];
    string colormap_fname = argv[2];
    string new_dir = argv[3];
    if(new_dir.back() == '/') new_dir = new_dir.substr(0,new_dir.size()-1);
    new_dir += '/';

    Mesh m;
    vector<Bem::real> phi,psi;
    import_ply(model_fname,m,phi,psi);

    // "/home/thomas/Documents/ma-thesis/code/first-tries/python-code/exterior-phi/binary.csv"
    ifstream cmap_input(colormap_fname);
    vector<Bem::real> pos;
    vector<vec3> cmap;
    string line;
    while(getline(cmap_input,line)) {
        stringstream linestr(line);
        Bem::real new_pos;
        vec3 new_col;
        linestr >> new_pos;
        pos.push_back(new_pos);
        linestr >> new_col.x;
        linestr >> new_col.y;
        linestr >> new_col.z;
        cmap.push_back(new_col);
    }

    cout << pos.size() << endl;

    Bem::real max_phi = phi[0];
    Bem::real min_phi = phi[0];

    if(argc == 6) {
        min_phi = atof(argv[4]);
        max_phi = atof(argv[5]);
    } else {
        for(Bem::real elm : phi) {
            max_phi = max(elm,max_phi);
            min_phi = min(elm,min_phi);
        }
    }

    vector<vec3> colors;
    for(Bem::real elm : phi) {
        Bem::real val = (elm - min_phi)/(max_phi-min_phi);
        colors.push_back(interpolate(pos,cmap,val));
    }
    

    cout << min_phi << ' ' << max_phi << endl;

    export_ply_colors(new_dir+model_fname.substr(model_fname.find_last_of('/')),m,colors);

}