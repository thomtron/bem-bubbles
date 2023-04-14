#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <string>
#include <fstream>
#include <cassert>

#include "../basic/Bem.hpp"
#include "Triplet.hpp"

/*
    General Note: Maybe it would be better to keep the Mesh class very slender, maybe make it even a struct;
    with public attributes and functions for initialisation only. all the other functions will then take a 
    Mesh reference or const Mesh reference as argument. I think this would improve adaptivity. And further 
    many other similar libraries work in the same way.
*/

namespace Bem {

class Tuplet {
public:
    Tuplet(size_t a_,size_t b_)
        :a(a_),b(b_) { 
            if(a_<b_){
                a = b_;
                b = a_;
            }
        }

    bool operator==(Tuplet const& other) {
        return (a==other.a and b==other.b); 
    }

    bool operator<(Tuplet const& other) { // "alphabetical sort"
        if(a<other.a) return true;
        if(a==other.a) return b<other.b;
        return false;
    }

    size_t get_a() const { return a; }
    size_t get_b() const { return b; }

private:
    size_t a,b;
};

struct Halfedge {
    Halfedge* twin;
    Halfedge* next;
    size_t vertex;
    size_t edge;
    size_t triangle;
};

struct HalfedgeMesh {
    std::vector<Halfedge*> verts;
    std::vector<Halfedge*> faces;
    std::vector<Halfedge*> edges;
    std::vector<Halfedge> halfs;
};


class Mesh {
public:
    using real = Bem::real;
    using vec3 = Bem::vec3;

    Mesh() {}

    virtual ~Mesh() {}

    void export_obj(std::string filename) const;
    void export_ply(std::string filename) const;
    void export_ply(std::string filename,std::vector<real> values,real min,real max) const;
    void export_ply_float(std::string filename,std::vector<real> values) const;
    void export_ply_double(std::string filename,std::vector<real> values) const;
    void export_ply_float_separat(std::string filename,std::vector<real> values) const;
    void export_ply_double_separat(std::string filename,std::vector<real> values) const;

    void import_ply(std::string filename,std::vector<real>& values);

    size_t N() const { return vertices.size(); }
    size_t M() const { return triangles.size(); }

    // index list of triangles that share the vertex i
    std::vector<std::vector<size_t>> generate_triangle_indices() const;
    std::vector<std::vector<size_t>> generate_2_ring() const;

    std::vector<std::vector<size_t>> generate_neighbours() const;

    std::vector<std::vector<size_t>> generate_edges() const;
    HalfedgeMesh generate_halfedges() const;
    std::vector<std::vector<size_t>> generate_triangle_edge_indices(std::vector<std::vector<size_t>> const& edges) const;
    
    void split_edges(real L_max);
    void collapse_edges(real L_min);
    void flip_edges(int state = 0); // 0-> valence and surface cost, 1-> valence cost, 2-> surface cost
    void relax_vertices();
    void project(Mesh const& other);
    std::vector<Bem::real> project_and_interpolate(Mesh const& other,std::vector<Bem::real> const& f);

    std::vector<vec3> generate_triangle_normals() const; // with vector product
    std::vector<vec3> generate_vertex_normals() const;   // weighted average of normals of triangles around vertex i
    std::vector<vec3> generate_vertex_normals_max() const;

    real solid_angle_at_vertex(size_t i) const; // solid angle of mesh w.r.t. vertex i
    real volume() const;
    real curvature() const;

    // creates a list of meshes, that are disjoint (share vertices but no triangles). 
    // Decomposition is based on indexing of triangles and vertices are stored for each (could be optimized therein)
    std::vector<Mesh> generate_mesh_partition(size_t num) const;

    const std::vector<vec3>& get_vertices() const {
        return vertices;
    }

    const std::vector<Triplet>& get_triangles() const {
        return triangles;
    }

    Triplet get_triangle(size_t i) const {
        return triangles[i];
    }

    // setters for vertices and triangles
    void set_vertices(std::vector<vec3> const& positions) {
        vertices = positions;
    }
    void set_triangles(std::vector<Triplet> const& trigs) {
        triangles = trigs;
    }


    void set_vertex_position(vec3 const& pos, size_t i) {
        vertices[i] = pos;
        return;
    }

    void add_vertex(vec3 const& pos) {
        vertices.push_back(pos);
    }

    void add_triangle(Triplet const& t) {
        
        triangles.push_back(t);
    }

    void add_triangle(size_t i,size_t j,size_t k) {
        
        add_triangle(Triplet(i,j,k));
    }

    void add_mesh(Mesh const& other,vec3 position = vec3(0.0,0.0,0.0));

    void clear() {
        vertices.clear();
        triangles.clear();
    }

protected:

    void trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index);

    // position in 3d space of the vertices of the mesh
    std::vector<vec3> vertices;
    // indices of the three vertices that build up the triangle
    std::vector<Triplet> triangles;  
};

template<typename T>
T parse_binary(std::ifstream& input) {
    std::vector<char> buf(sizeof(T));
    input.read(buf.data(),buf.size());
    T* res = reinterpret_cast<T*>(buf.data());
    return *res;
}

template<typename T>
T parse_binary(std::ifstream& input,size_t buffersize) {
    std::vector<char> buf(sizeof(T));
    //assert(buffersize>=buf.size());
    input.read(buf.data(),buffersize);
    T* res = reinterpret_cast<T*>(buf.data());
    return *res;
}

} // namespace Bem


#endif // MESH_HPP
