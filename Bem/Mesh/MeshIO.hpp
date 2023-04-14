#ifndef MESHIO_HPP
#define MESHIO_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../basic/Bem.hpp"
#include "Mesh.hpp"

namespace Bem {

void export_obj                 (std::string filename, Mesh const& mesh);
void export_ply                 (std::string filename, Mesh const& mesh);
void export_ply                 (std::string filename, Mesh const& mesh, std::vector<real> const& values, real min, real max);
void export_ply_colors          (std::string filename, Mesh const& mesh, std::vector<vec3> const& colors);
void export_ply_float           (std::string filename, Mesh const& mesh, std::vector<real> const& values);
void export_ply_double          (std::string filename, Mesh const& mesh, std::vector<real> const& values);
void export_ply_float           (std::string filename, Mesh const& mesh, std::vector<real> const& phi, std::vector<real> const& psi);
void export_ply_double          (std::string filename, Mesh const& mesh, std::vector<real> const& phi, std::vector<real> const& psi);
void export_ply_float_separat   (std::string filename, Mesh const& mesh, std::vector<real> const& values);
void export_ply_double_separat  (std::string filename, Mesh const& mesh, std::vector<real> const& values);

void import_ply                 (std::string filename, Mesh& mesh, std::vector<real>& values);
void import_ply                 (std::string filename, Mesh& mesh, std::vector<real>& phi, std::vector<real>& psi);
void import_ply                 (std::string filename, Mesh& mesh);

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