#ifndef MESHMANIP_HPP
#define MESHMANIP_HPP

#include "Mesh.hpp"
#include "HalfedgeMesh.hpp"


// following code from: https://stackoverflow.com/questions/44929500/how-to-get-a-stdlisttiterator-from-an-element-of-that-list
// probably not good practice!

//This is essentially what you are looking for:
template<typename T>
typename std::list<T>::iterator pointerToIter (T* myPointer) {
    //Calculates the distance in bytes from an iterator itself
    //to the actual type that is stored at the position the
    //iterator is pointing to.
    size_t iterOffset = (size_t)&(*((std::list<void*>::iterator)nullptr));
    //Subtract the offset from the passed pointer and make an
    //iterator out of it
    typename std::list<T>::iterator iter;
    *(intptr_t*)&iter = (intptr_t)myPointer - iterOffset;
    //You are done
    return iter;
}

namespace Bem {

// This file contains a collection of functions that are useful for manipulating a triangle mesh.
// split-, collapse- and flip_edges are connectivity-changing remeshing functions. split- and
// collapse_edges further add resp. remove elements from the mesh. relax_vertices is a smoothing
// function that only displaces vertices. the second group of functions are needed for projecting
// a mesh along its normals on another mesh and interpolating vertex data from the other mesh.

void split_edges    (HalfedgeMesh& mesh, real L_max);
void split_edges    (HalfedgeMesh& mesh, std::vector<real>& max_edgelenght);
void collapse_edges (HalfedgeMesh& mesh, real L_min);
void collapse_edges (HalfedgeMesh& mesh, std::vector<real>& curvature, real multiplicator);
void flip_edges     (HalfedgeMesh& mesh, size_t state = 0);
void relax_vertices (Mesh& mesh);
void relax_vertices (HalfedgeMesh& mesh);

// todo: let trace_mesh return a bool too.
void trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index);
bool trace_mesh_positive(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index);
void project(Mesh& mesh, Mesh const& other);
void project_from_origin(std::vector<vec3>& normals, Mesh const& other, real const& dist_to_wall)
void project_and_interpolate(Mesh& mesh, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);
void project_and_interpolate(Mesh& mesh, std::vector<vec3> const& vertex_normals, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);

Mesh l2smooth(Mesh mesh);
Mesh l2smooth(Mesh mesh, std::vector<size_t> const& vert_inds);
Mesh l2smooth(Mesh mesh, std::vector<real>& pot);
Mesh l2smooth(Mesh mesh, std::vector<real>& pot, std::vector<size_t> const& vert_inds);

void check_validity(HalfedgeMesh const& mesh);

} // namespace Bem

#endif // MESHMANIP_HPP