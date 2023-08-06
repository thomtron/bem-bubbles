#ifndef HALFEDGEMESH_HPP
#define HALFEDGEMESH_HPP

#include <iostream>
#include <list>

#include "Mesh.hpp"

namespace Bem {

// The HalfedgeMesh is a structure describing a triangle mesh by means of so called Halfedges.
// Each edge has two Halfedges that point in oposite directions. A halfedge points always from one
// vertex to the other and it has a pointer to its twin, to the next halfedge (rooting at the vertex
// it points to) of the same triangle, to a vertex (where its root is), to its edge and its triangle.
// To form a complete structure, the vertices, edges and triangles each have to have a pointer to one 
// Halfedge they're connected with too. Although it isn't important which exact Halfedge they point to.
// This structure is quite a bit more complicated, but it allows to jump easily from one neighbour of 
// vertices/edges/triangles to another without the need to ever loop over all elements. Most applications
// of this structure are presented in MeshManip.hpp/.cpp. In HalfedgeMesh.cpp we define the functions 
// needed to translate between the Mesh and HalfedgeMesh representations.

struct Vertex;

struct Halfedge {
    Halfedge* twin;
    Halfedge* next;
    //size_t vert;
    //size_t edge;
    //size_t trig;
    Vertex* vert;
    Halfedge** edge;
    Halfedge** trig;

};

struct Vertex {
    Halfedge* half;
    vec3 pos;
    size_t index;
};

struct HalfedgeMesh {
    std::list<Vertex>    verts;
    std::list<Halfedge*> trigs;
    std::list<Halfedge*> edges;
    std::list<Halfedge*> bounds;

    ~HalfedgeMesh();

    void release();

    /*
    template<typename T>
    size_t get_index(std::vector<T> const& vec, const T* elm) const {
        size_t index = elm - &vec[0]; // internet says: ()/sizeof(T), but that didn't give correct results...
        //std::cout << elm << " + " << &vec[0] << " + " << index << std::endl;
        if(index >= vec.size()) throw(std::out_of_range("Bad Vertex pointer: out of range of Vertex vector."));
        return index;
    }*/
    
    HalfedgeMesh(Mesh const& other);
    HalfedgeMesh(HalfedgeMesh&&) = default;
    HalfedgeMesh() = default;
    HalfedgeMesh(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh&&) = default;
    
    void copy(HalfedgeMesh const& other);

    void update_vertex_indices();
};

HalfedgeMesh generate_halfedges(Mesh const& mesh);
void generate_halfedges(HalfedgeMesh& result, Mesh const& mesh);
Mesh generate_mesh(HalfedgeMesh const& mesh);
void generate_mesh(Mesh& result, HalfedgeMesh const& mesh);

} // namespace Bem

#endif // HALFEDGEMESH_HPP