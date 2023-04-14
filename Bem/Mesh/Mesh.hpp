#ifndef MESH_NEW_HPP
#define MESH_NEW_HPP

#include <vector>
#include <string>
#include <cassert>

#include "../basic/Bem.hpp"
#include "Triplet.hpp"

namespace Bem {


struct Mesh {
    // position in 3d space of the vertices of the mesh
    std::vector<vec3> verts;
    // indices of the three vertices that build up the triangle
    std::vector<Triplet> trigs;  

    void add(Mesh other,vec3 const& position = vec3());
    void scale(real s);
    void rotate(vec3 v);
    void clear();
};

std::vector<std::vector<size_t>> generate_triangle_indices(Mesh const& mesh);
std::vector<std::vector<size_t>> generate_2_ring(Mesh const& mesh);
std::vector<std::vector<size_t>> generate_2_ring(Mesh const& mesh, std::vector<std::vector<size_t>> const& neighbours);
std::vector<std::vector<size_t>> generate_neighbours(Mesh const& mesh);

std::vector<vec3> generate_triangle_normals  (Mesh const& mesh); // with vector product
std::vector<vec3> generate_vertex_normals    (Mesh const& mesh, std::vector<std::vector<size_t>> const& triangle_indices);
std::vector<vec3> generate_vertex_normals    (Mesh const& mesh);

// solid angle of mesh w.r.t. vertex i
real solid_angle_at_vertex(Mesh const& mesh, std::vector<std::vector<size_t>> const& triangle_indices, std::vector<vec3>& triangle_normals, size_t i); 
real solid_angle_at_vertex(Mesh const& mesh, size_t i);
real volume(Mesh const& mesh); // computes total volume of the mesh
vec3 centerofmass(Mesh const& mesh); // computes center of mass for constant mass density
void to_centerofmass(Mesh& mesh);

std::vector<real> max_curvature(Mesh const& mesh);
void curvatures(Mesh const& mesh, std::vector<real>& kappa, std::vector<real>& gamma);

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
