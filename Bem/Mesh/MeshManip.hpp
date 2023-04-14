#ifndef MESHMANIP_HPP
#define MESHMANIP_HPP

#include "Mesh.hpp"
#include "HalfedgeMesh.hpp"

namespace Bem {

class CubicInterpolator {
public:
    CubicInterpolator(vec3 a,vec3 b,vec3 c,vec3 na,vec3 nb, vec3 nc)
        :p300(a),p030(b),p003(c) {
            p210 = project_twothrids(a,b,na);
            p120 = project_twothrids(b,a,nb);
            p021 = project_twothrids(b,c,nb);
            p012 = project_twothrids(c,b,nc);
            p102 = project_twothrids(c,a,nc);
            p201 = project_twothrids(a,c,na);
            vec3 E = (p210 + p120 + p021 + p012 + p102 + p201)*(1.0/6.0);
            vec3 V = (a+b+c)*(1.0/3.0);
            p111 = E + (E-V)*0.5;  
        }

    vec3 interpolate(real u,real v) const {
        real w = 1.0 - u - v;
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        
        return u*u2*p300 + v*v2*p030 + w*w2*p003 + 3*u2*v*p210 + 3*u2*w*p201 + 3*u*v2*p120 + 3*v2*w*p021 + 3*v*w2*p012 + 3*u*w2*p102 + 6*u*v*w*p111;
    }

    vec3 get_surface_vector(real u,real v) const {
        real w = 1.0 - u - v;
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        vec3 dfdu = 3.0*(u2*p300 - w2*p003 + v2*(p120-p021)) + (6.0*u*w - 3.0*u2)*p201 + (3.0*w2 - 6.0*u*w)*p102 + 6.0*((v*w - u*v)*p111 + u*v*p210 - v*w*p012);
        vec3 dfdv = 3.0*(v2*p030 - w2*p003 + u2*(p210-p201)) + (6.0*v*w - 3.0*v2)*p021 + (3.0*w2 - 6.0*v*w)*p012 + 6.0*((u*w - u*v)*p111 + u*v*p120 - u*w*p102);

        return dfdu.vec(dfdv);
    }

    vec3 get_normal(real u,real v) const {
        vec3 normal(get_surface_vector(u,v));
        normal.normalize();
        return normal;
    }

private:

    vec3 project_twothrids(vec3 const& p1,vec3 const& p2,vec3 const& n1) const {
        real w = (p1-p2).dot(n1);
        return (2.0*p1 + p2 + w*n1)*(1.0/3.0);
    }

    vec3 p300,p030,p003,p210,p120,p021,p012,p102,p201,p111;
};

void split_edges    (HalfedgeMesh& mesh, real L_max);
void split_edges    (HalfedgeMesh& mesh, std::vector<real>& curvature, real multiplicator);
void collapse_edges (HalfedgeMesh& mesh, real L_min);
void collapse_edges (HalfedgeMesh& mesh, std::vector<real>& curvature, real multiplicator);
void flip_edges     (HalfedgeMesh& mesh, size_t state = 0);
void relax_vertices (Mesh& mesh);
void relax_vertices (HalfedgeMesh& mesh);

void trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index);
void project   (Mesh& mesh, Mesh const& other);
void project_and_interpolate(Mesh& mesh, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);
void project_and_interpolate(Mesh& mesh, std::vector<vec3> const& vertex_normals, std::vector<real>& f_res, Mesh const& other, std::vector<real> const& f);

} // namespace Bem

#endif // MESHMANIP_HPP