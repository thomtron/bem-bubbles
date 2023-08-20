#ifndef MESH_NEW_HPP
#define MESH_NEW_HPP

#include <vector>
#include <string>
#include <cassert>

#include "../basic/Bem.hpp"

namespace Bem {

// the mesh class consists of a vector of vertices and a vector triangle indices.
// additionaly there are some very basic class methods.

struct Mesh {
    // position in 3d space of the vertices of the mesh
    std::vector<vec3> verts;
    // indices of the three vertices that build up the triangle
    std::vector<Triplet> trigs;  

    void add(Mesh other,vec3 const& position = vec3());
    void scale(real s);
    void rotate(vec3 v);
    void translate(vec3 v);
    void clear();

    bool check_validity() const;
};

// generate_triangle_indices returns for each vertex a vector of indices 
// of the triangles touching this vertex. The outer vector's length is thus
// the same as mesh.verts.size()
std::vector<std::vector<size_t>> generate_triangle_indices(Mesh const& mesh);
// generate_neighbours returns for each vertex a list of its direct neighbours (connected with an edge)
std::vector<std::vector<size_t>> generate_neighbours(Mesh const& mesh);
// generate_2_ring returns for each vertex a list of the neighbours of its neighbours, including the vertex itself
std::vector<std::vector<size_t>> generate_2_ring(Mesh const& mesh);
std::vector<std::vector<size_t>> generate_2_ring(Mesh const& mesh, std::vector<std::vector<size_t>> const& neighbours);

std::vector<vec3> generate_triangle_normals  (Mesh const& mesh); // with vector product
// vertex normals according to Max_1999
std::vector<vec3> generate_vertex_normals    (Mesh const& mesh, std::vector<std::vector<size_t>> const& triangle_indices);
std::vector<vec3> generate_vertex_normals    (Mesh const& mesh);

// solid angle of mesh w.r.t. vertex i
real solid_angle_at_vertex(Mesh const& mesh, std::vector<std::vector<size_t>> const& triangle_indices, std::vector<vec3>& triangle_normals, size_t i); 
real solid_angle_at_vertex(Mesh const& mesh, size_t i);

real volume(Mesh const& mesh);          // computes total volume of the mesh

vec3 centerofmass(Mesh const& mesh);    // computes center of mass for constant mass density
void to_centerofmass(Mesh& mesh);       // translates the mesh such that its center of mass is at (0,0,0)

// functions that compute the curvature of the mesh
std::vector<real> max_curvature(Mesh const& mesh); // maximum curvature
void curvatures(Mesh const& mesh, std::vector<real>& kappa, std::vector<real>& gamma); // mean and gaussian curvature


// the following functions are useful for the simulation of multiple bubbles that are
// for example initialized with different potentials or whose volume has to be computed
// separately for time evolution with nonzero gas pressure

Mesh join_meshes(std::vector<Mesh> const& list);
std::vector<Mesh> split_by_loose_parts(Mesh const& mesh);
std::vector<Mesh> split_by_loose_parts(Mesh const& mesh, std::vector<std::vector<size_t>>& vert_perm);

// VD stands for Vertex Data
template<typename T>
std::vector<T> expand_VD_to_joined(std::vector<Mesh> const& group, std::vector<std::vector<T>> const& separated) {
    assert(group.size() == separated.size());
    for(size_t i(0);i<group.size();++i)
        assert(group[i].verts.size() == separated[i].size());

    std::vector<T> result;
    for(size_t i(0);i<group.size();++i)
        for(size_t j(0);j<group[i].verts.size();++j)
            result.push_back(separated[i][j]);

    return result;
}

template<typename T>
std::vector<std::vector<T>> collapse_VD_to_separated(std::vector<Mesh> const& group, std::vector<T> const& joined) {
    size_t total_size(0);
    for(Mesh const& elm : group)
        total_size += elm.verts.size();
    assert(total_size == joined.size());

    std::vector<std::vector<T>> result;
    size_t offset(0);
    for(Mesh const& elm : group) {
        std::vector<T> part;
        for(size_t i(0);i<elm.verts.size();++i) {
            part.push_back(joined[i + offset]);
        }
        offset += elm.verts.size();
        result.push_back(part);
    }
    return result;
}


} // namespace Bem


#endif // MESH_NEW_HPP
