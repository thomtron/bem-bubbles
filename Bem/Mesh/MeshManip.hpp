#ifndef MESHMANIP_HPP
#define MESHMANIP_HPP

#include "Mesh.hpp"
#include "HalfedgeMesh.hpp"

namespace Bem {

// This file contains a collection of functions that are useful for manipulating a triangle mesh.
// split-, collapse- and flip_edges are connectivity-changing remeshing functions. split- and
// collapse_edges further add resp. remove elements from the mesh. relax_vertices is a smoothing
// function that only displaces vertices. the second group of functions are needed for projecting
// a mesh along its normals on another mesh and interpolating vertex data from the other mesh.

void split_edges    (HalfedgeMesh& mesh, real L_max);
void split_edges    (HalfedgeMesh& mesh, std::vector<real>& curvature, real multiplicator);
void collapse_edges (HalfedgeMesh& mesh, real L_min);
void collapse_edges (HalfedgeMesh& mesh, std::vector<real>& curvature, real multiplicator);
void flip_edges     (HalfedgeMesh& mesh, size_t state = 0);
void relax_vertices (Mesh& mesh);
void relax_vertices (HalfedgeMesh& mesh);

bool trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index);
void project   (Mesh& mesh, Mesh const& other);
void project_and_interpolate(Mesh& mesh, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);
void project_and_interpolate(Mesh& mesh, std::vector<vec3> const& vertex_normals, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);
void project_and_interpolate(Mesh& mesh, std::vector<vec3> const& vertex_normals, std::vector<real>& f_res, std::vector<real>& f_2_res, Mesh const& other, std::vector<real> const& f, std::vector<real> const& f_2);

Mesh l2smooth(Mesh mesh);
Mesh l2smooth(Mesh mesh, std::vector<size_t> const& vert_inds);
Mesh l2smooth(Mesh mesh, std::vector<real>& pot);
Mesh l2smooth(Mesh mesh, std::vector<real>& pot, std::vector<size_t> const& vert_inds);

} // namespace Bem

#endif // MESHMANIP_HPP