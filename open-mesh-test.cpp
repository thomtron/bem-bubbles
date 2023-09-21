#include <iostream>
#include <vector>

#include "Bem/basic/Bem.hpp"
#include "Bem/Mesh/Mesh.hpp"
#include "Bem/Mesh/MeshIO.hpp"
#include "Bem/Mesh/MeshManip.hpp"

using namespace std;
using namespace Bem;

int main() {
    Mesh M;
    import_ply("semi-sphere.ply",M);

    
    HalfedgeMesh half; //(generate_halfedges(M));

    generate_halfedges(half,M);
    if(half.check_validity()) cout << "Hurray!" << endl;
    
    Bem::real L(0.2);

    vector<Bem::real> curvs(half.edges.size(),1.0);


    split_edges(half,curvs,L*4.0/3.0);
    collapse_edges(half,curvs,L*4.0/5.0);
    split_edges(half,curvs,L*4.0/3.0);
    collapse_edges(half,curvs,L*4.0/5.0);
    split_edges(half,curvs,L*4.0/3.0);
    collapse_edges(half,curvs,L*4.0/5.0);
    flip_edges(half,1);
    flip_edges(half,1);

    relax_vertices(half);
    export_ply("split_edges.ply",generate_mesh(half));
    /*
    split_edges(half,L*4.0/3.0);
    collapse_edges(half,L*4.0/5.0);
    split_edges(half,L*4.0/3.0);
    collapse_edges(half,L*4.0/5.0);
    split_edges(half,L*4.0/3.0);
    collapse_edges(half,L*4.0/5.0);
    split_edges(half,L*4.0/3.0);
    collapse_edges(half,L*4.0/5.0);
    split_edges(half,L*4.0/3.0);
    collapse_edges(half,L*4.0/5.0);
    flip_edges(half,1);
    flip_edges(half,1);

    relax_vertices(half);*/


    //flip_edges(half,1);
    //flip_edges(half,1);
    //relax_vertices(half);
    //collapse_edges(half,L*4.0/5.0);
    
    if(half.check_validity()) cout << "Hurray!" << endl;
    export_ply("split_edges.ply",generate_mesh(half));



    return 0;
}