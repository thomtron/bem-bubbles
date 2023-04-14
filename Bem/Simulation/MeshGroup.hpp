#ifndef MESHGROUP_HPP
#define MESHGROUP_HPP

#include <iostream>
#include <vector>
#include "../Mesh/Mesh.hpp"


namespace Bem {

class MeshGroup {
public:
    MeshGroup(std::vector<Mesh> const& meshes);
    MeshGroup() = default;
    void add(std::vector<Mesh> const& meshes);
    void add(Mesh const& mesh);

    Mesh const& get_joined() {
        update_joined();
        return joined;
    }

    // returns for each vertex the volume of its corresponding mesh
    std::vector<real> volumes() const;
    std::vector<real> volume_per_vertex() const;
    std::vector<real> expand_to_vertex_data(std::vector<real> const& bubble_data) const;
    std::vector<std::vector<real>> collapse_to_bubble_data(std::vector<real> const& vertex_data) const;
    
    size_t num_verts() const;
    size_t num_trigs() const;

    void set_vertices(std::vector<vec3> verts) {
        assert(num_verts() == verts.size());
        size_t offset = 0;
        for(Mesh& mesh : separat) {
            for(size_t i(0);i<mesh.verts.size();++i)
                mesh.verts[i] = verts[offset+i];
            offset = mesh.verts.size();
        }
    }

private:
    // update functions - They only make a change if mod_joined/mod_separat are true
    // They are only called in functions where the respective attribute is needed
    void update_joined();
    void update_separat();

    Mesh joined;
    std::vector<Mesh> separat;
    bool recent; // if joined is most recent version, recent = 1, otherwise recent = 0
};



} // namespace Bem

#endif // MESHGROUP_HPP