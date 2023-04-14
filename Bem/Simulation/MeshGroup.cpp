#include "MeshGroup.hpp"

using namespace std;

namespace Bem {

MeshGroup::MeshGroup(vector<Mesh> const& meshes) {
    add(meshes);
    update_joined();
}

void MeshGroup::add(Mesh const& mesh) {
    separat.push_back(mesh);
}

void MeshGroup::add(vector<Mesh> const& meshes) {
    for(Mesh const& mesh : meshes) separat.push_back(mesh);
}

vector<real> MeshGroup::volumes() const {
    vector<real> result;
    for(Mesh const& mesh : separat) {
        result.push_back(volume(mesh));
    }
    return result;
}

vector<real> MeshGroup::volume_per_vertex() const {
    return expand_to_vertex_data(volumes());
}

vector<real> MeshGroup::expand_to_vertex_data(vector<real> const& bubble_data) const {
    assert(bubble_data.size() == separat.size());
    vector<real> result;
    for(size_t i(0);i<separat.size();++i)
        for(size_t j(0);j<separat[i].verts.size();++j)
            result.push_back(bubble_data[i]);

    return result;
}

vector<vector<real>> MeshGroup::collapse_to_bubble_data(vector<real> const& vertex_data) const {
    assert(num_verts() == vertex_data.size());
    vector<vector<real>> result;
    size_t offset(0);
    for(Mesh const& mesh : separat) {
        vector<real> bubble_dat;
        for(size_t i(0);i<mesh.verts.size();++i) {
            bubble_dat.push_back(vertex_data[i + offset]);
        }
        result.push_back(bubble_dat);
    }
    return result;
}

size_t MeshGroup::num_verts() const {
    if(recent == 0) {
        size_t res(0);
        for(Mesh const& mesh : separat)
            res += mesh.verts.size();
        return res;
    } else {
        return joined.verts.size();
    }
}
size_t MeshGroup::num_trigs() const {
    if(recent == 0) {
        size_t res(0);
        for(Mesh const& mesh : separat)
            res += mesh.trigs.size();
        return res;
    } else {
        return joined.trigs.size();
    }
}

void MeshGroup::update_joined() {
    if(recent == 0) {
        recent = 1;
        joined.clear();
        for(Mesh const& mesh : separat) joined.add(mesh);
    }
}

void MeshGroup::update_separat() {
    if(recent == 1) {
        recent = 0;

        size_t n(0);
        for(Mesh const& mesh : separat) n += mesh.verts.size();
        assert(n == joined.verts.size());

        // update verts of separat
        size_t offset(0);
        for(Mesh& mesh : separat) {
            for(size_t i(0);i<mesh.verts.size();++i) {
                mesh.verts[i] = joined.verts[offset+i];
            }
            offset += mesh.verts.size();
            // Note, the Topology of the meshes must not be changed in joined, 
            // so we can leave the trigs as they are!
        }
    }
}

} // namespace Bem