#include <vector>
#include <algorithm> // for sort()
#include "Icosphere.hpp"

using namespace std;

using namespace Bem;

// tuplet is a small helper struct for representing edges. It is only used by the
// create_icosphere function. May be replaced by standard library type?
struct tuplet {
public:
    size_t a,b;
    tuplet(size_t a_,size_t b_)
        :a(a_),b(b_) { 
            if(a_>b_){
                a = b_;
                b = a_;
            }
        }

    bool operator==(tuplet const& other) {
        return (a==other.a and b==other.b) or (a==other.b and b==other.a); 
    }

    bool operator<(tuplet const& other) { // "alphabetical sort"
        if(a<other.a) return true;
        if(a==other.a) return b<other.b;
        return false;
    }

    // make a hash: i = 2**a * 3**b or something like an id-generator (static object)
};

void Icosphere::create_icosphere(unsigned int order) {

    const real phi(0.5*(1.0+sqrt(5.0)));

    vertices.clear();
    triangles.clear();

    // Quelle angeben!!
    vertices.push_back(vec3(0.0, 1.0, phi));
    vertices.push_back(vec3(0.0,-1.0, phi)); 
    vertices.push_back(vec3(1.0, phi, 0.0));
    vertices.push_back(vec3(-1.0, phi, 0.0)); 
    vertices.push_back(vec3(phi, 0.0, 1.0)); 
    vertices.push_back(vec3(-phi, 0.0, 1.0));

    // rescale vertices (so that they're normalized) and add the negative vectors
    for(size_t i(0);i<6;++i){
        vertices[i] *= 1.0/sqrt(phi*phi + 1.0);
        vertices.push_back(-vertices[i]);
    }

    // vertices of basis icosahedron are now ready. Next come the faces:
    triangles.push_back(Triplet(0,5,1));
    triangles.push_back(Triplet(0,3,5));
    triangles.push_back(Triplet(0,2,3));
    triangles.push_back(Triplet(0,4,2));
    triangles.push_back(Triplet(0,1,4));
    triangles.push_back(Triplet(1,5,8));
    triangles.push_back(Triplet(5,3,10));
    triangles.push_back(Triplet(3,2,7));
    triangles.push_back(Triplet(2,4,11));
    triangles.push_back(Triplet(4,1,9));
    triangles.push_back(Triplet(7,11,6));
    triangles.push_back(Triplet(11,9,6));
    triangles.push_back(Triplet(9,8,6));
    triangles.push_back(Triplet(8,10,6));
    triangles.push_back(Triplet(10,7,6));
    triangles.push_back(Triplet(2,11,7));
    triangles.push_back(Triplet(4,9,11));
    triangles.push_back(Triplet(1,8,9));
    triangles.push_back(Triplet(5,10,8));
    triangles.push_back(Triplet(3,7,10));

    
    for(unsigned int k(0); k<order; ++k){

        vector<tuplet> edges;
        vector<Triplet> edge_indices;
        size_t ind(0);
        for(Triplet const& tri : triangles){
            edges.push_back(tuplet(tri.a,tri.b));
            edges.push_back(tuplet(tri.b,tri.c));
            edges.push_back(tuplet(tri.c,tri.a));
            edge_indices.push_back(Triplet(ind,ind+1,ind+2));
            ind += 3;
        }

        // each edge has been counted twice -> we'll sort the edges and extract each second
        
        // create an index list of the size of edges
        vector<size_t> indices;
        for(size_t i(0);i<edges.size();++i){
            indices.push_back(i);
        }
        // sort the indices (w.r.t. the alphabetical order of the edges -> remove duplicates)
        sort(indices.begin(),indices.end(),[&edges](size_t a,size_t b) { return edges[a] < edges[b]; });

        vector<tuplet> new_edges;
        for(size_t i(0);i<edges.size();i+=2){
            size_t ind = indices[i];
            new_edges.push_back(edges[ind]);

            // update the edge indices (of each triangle)
            edge_indices[ind/3][ind%3] = i/2;
            ind = indices[i+1];
            edge_indices[ind/3][ind%3] = i/2; // new index corresponds to index in new_edges vector
            // now we have for each triangle, the indices of the corresponding edges.
        }
        edges.clear();
        edges = new_edges;
        new_edges.clear();



        // test one refinement:

        size_t first_index(vertices.size());
        for(tuplet const& t : edges){
            vertices.push_back(vertices[t.a] + vertices[t.b]);
            vertices.back().normalize();
        }

        vector<Triplet> new_triangles;
        for(size_t i(0);i<triangles.size();++i){
            Triplet t(triangles[i]);
            Triplet ind(edge_indices[i]);
            ind.a += first_index;
            ind.b += first_index;
            ind.c += first_index;
            new_triangles.push_back(Triplet(t.a,ind.a,ind.c));
            new_triangles.push_back(Triplet(ind.a,t.b,ind.b));
            new_triangles.push_back(Triplet(ind.c,ind.b,t.c));
            new_triangles.push_back(Triplet(ind.a,ind.b,ind.c));
        }
        triangles.clear();
        triangles = new_triangles;
        new_triangles.clear();

    }

}