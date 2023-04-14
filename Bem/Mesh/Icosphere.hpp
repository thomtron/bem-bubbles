#ifndef ICOSPHERE_HPP
#define ICOSPHERE_HPP

#include "Mesh.hpp"

namespace Bem {

class Icosphere : public Mesh {
public:
    Icosphere(unsigned int order) {
        create_icosphere(order);

        // experimental:
        
        /*
        size_t n(N());
        for(size_t i(0);i<n;++i){
            vertices.push_back(vec3(2.3,0.0,0.0)+vertices[i]); // mirrored sphere
        }
        size_t m(M());
        for(size_t i(0);i<m;++i){
            triangles.push_back(Triplet(triangles[i].a+n,triangles[i].b+n,triangles[i].c+n));
        }
        */
    }
    virtual ~Icosphere() {}
private:
    void create_icosphere(unsigned int order);
};

} // namespace Bem

#endif // ICOSPHERE_HPP