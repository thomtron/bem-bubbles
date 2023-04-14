#ifndef HALFEDGEMESH_HPP
#define HALFEDGEMESH_HPP

#include <iostream>
#include <vector>

#include "Mesh.hpp"

namespace Bem {

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

    ~HalfedgeMesh();
    
    HalfedgeMesh(Mesh const& other);
    HalfedgeMesh(HalfedgeMesh&&) = default;
    HalfedgeMesh() = default;
    HalfedgeMesh(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh const& other);
    HalfedgeMesh& operator=(HalfedgeMesh&&) = default;
    
    void copy(HalfedgeMesh const& other);
};

HalfedgeMesh generate_halfedges(Mesh const& mesh);
Mesh generate_mesh(HalfedgeMesh const& mesh);

} // namespace Bem

#endif // HALFEDGEMESH_HPP