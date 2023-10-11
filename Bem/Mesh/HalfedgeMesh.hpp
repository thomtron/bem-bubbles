#ifndef HALFEDGEMESH_HPP
#define HALFEDGEMESH_HPP

#include <iostream>
#include <vector>

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

struct Halfedge {
    Halfedge* twin;
    Halfedge* next;
    size_t vert;
    size_t edge;
    size_t trig;
};

struct HalfedgeMesh {
    std::vector<vec3>      vpos;
    std::vector<Halfedge*> verts;
    std::vector<Halfedge*> trigs;
    std::vector<Halfedge*> edges;
    std::vector<Halfedge*> bounds;

    ~HalfedgeMesh();

    void clear();
    bool check_validity() const;
    
    HalfedgeMesh(Mesh const& other);
    HalfedgeMesh(HalfedgeMesh&&) = default;
    HalfedgeMesh() = default;
    HalfedgeMesh(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh&&) = default; // check this!
    
    void copy(HalfedgeMesh const& other);

    static const size_t npos;
};

HalfedgeMesh generate_halfedges(Mesh const& mesh);
void generate_halfedges(HalfedgeMesh& result, Mesh const& mesh);
Mesh generate_mesh(HalfedgeMesh const& mesh);
void generate_mesh(Mesh& result, HalfedgeMesh const& mesh);

std::vector<vec3> generate_vertex_normals(HalfedgeMesh const& mesh);

} // namespace Bem

#endif // HALFEDGEMESH_HPP