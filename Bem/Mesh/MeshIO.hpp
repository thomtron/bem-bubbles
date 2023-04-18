#ifndef MESHIO_HPP
#define MESHIO_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../basic/Bem.hpp"
#include "Mesh.hpp"

namespace Bem {

// functions for the export of 3d meshes

void export_obj                 (std::string filename, Mesh const& mesh);                                                               // export as obj
void export_ply                 (std::string filename, Mesh const& mesh);                                                               // export as ply
void export_ply                 (std::string filename, Mesh const& mesh, std::vector<real> const& values, real min, real max);          // ply with vertex colors fading from blue to red between min and max (useful for quick visualisation)
void export_ply_colors          (std::string filename, Mesh const& mesh, std::vector<vec3> const& colors);                              // ply with vertex colors
void export_ply_float           (std::string filename, Mesh const& mesh, std::vector<real> const& values);                              // ply with a floating point value per vertex
void export_ply_double          (std::string filename, Mesh const& mesh, std::vector<real> const& values);                              // same with double precision                
void export_ply_float           (std::string filename, Mesh const& mesh, std::vector<real> const& phi, std::vector<real> const& psi);   // ply with two floating point values phi and psi per vertex
void export_ply_double          (std::string filename, Mesh const& mesh, std::vector<real> const& phi, std::vector<real> const& psi);   // same with double precision
void export_ply_float_separat   (std::string filename, Mesh const& mesh, std::vector<real> const& values);                              // ply with three floating point values per triangle
void export_ply_double_separat  (std::string filename, Mesh const& mesh, std::vector<real> const& values);                              // same with double precision

// functions for the import of 3d meshes

void import_ply                 (std::string filename, Mesh& mesh);                                                                     // imports ply
void import_ply                 (std::string filename, Mesh& mesh, std::vector<real>& values);                                          // imports ply with one additional vertex value
void import_ply                 (std::string filename, Mesh& mesh, std::vector<real>& phi, std::vector<real>& psi);                     // imports ply with two additional vertex values


// The following two templates are used to parse the binary data part of the ply files

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

#endif // MESHIO_HPP