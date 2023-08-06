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

    HalfedgeMesh half(generate_halfedges(M));

    Halfedge** f = half.trigs.back()->trig;
    cout << half.get_index(half.trigs,f) << endl << endl;

    Halfedge** h = half.edges[134]->edge;
    cout << half.get_index(half.edges,h) << endl << endl;

    

    Bem::real L(0.1);

    check_validity(half);
    split_edges(half,L*4.0/3.0);
    //flip_edges(half,1);
    //flip_edges(half,1);
    //relax_vertices(half);
    //collapse_edges(half,L*4.0/5.0);

    export_ply("split_edges.ply",generate_mesh(half));



    return 0;
}