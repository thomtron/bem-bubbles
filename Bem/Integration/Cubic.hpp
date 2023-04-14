#ifndef CUBIC_HPP
#define CUBIC_HPP

#include "../basic/Bem.hpp"

namespace Bem {


class Cubic {
public:
    Cubic(vec3 a,vec3 b,vec3 c,vec3 na,vec3 nb, vec3 nc)
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
        //if(u<0.0 or v<0.0 or w<0.0) std::cout << "error!";
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        
        return u*u2*p300 + v*v2*p030 + w*w2*p003 + 3*u2*v*p210 + 3*u2*w*p201 + 3*u*v2*p120 + 3*v2*w*p021 + 3*v*w2*p012 + 3*u*w2*p102 + 6*u*v*w*p111;
    }

    vec3 get_dudx(real u,real v) const {
        real w = 1.0 - u - v;
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        return 3.0*(u2*p300 - w2*p003 + v2*(p120-p021)) + (6.0*u*w - 3.0*u2)*p201 + (3.0*w2 - 6.0*u*w)*p102 + 6.0*((v*w - u*v)*p111 + u*v*p210 - v*w*p012);
    }
    vec3 get_dvdx(real u,real v) const {
        real w = 1.0 - u - v;
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        return 3.0*(v2*p030 - w2*p003 + u2*(p210-p201)) + (6.0*v*w - 3.0*v2)*p021 + (3.0*w2 - 6.0*v*w)*p012 + 6.0*((u*w - u*v)*p111 + u*v*p120 - u*w*p102);
    }

    vec3 get_surface_vector(real u,real v) const {
        real w = 1.0 - u - v;
        real u2 = u*u;
        real v2 = v*v;
        real w2 = w*w;
        vec3 dudx = 3.0*(u2*p300 - w2*p003 + v2*(p120-p021)) + (6.0*u*w - 3.0*u2)*p201 + (3.0*w2 - 6.0*u*w)*p102 + 6.0*((v*w - u*v)*p111 + u*v*p210 - v*w*p012);
        vec3 dvdx = 3.0*(v2*p030 - w2*p003 + u2*(p210-p201)) + (6.0*v*w - 3.0*v2)*p021 + (3.0*w2 - 6.0*v*w)*p012 + 6.0*((u*w - u*v)*p111 + u*v*p120 - u*w*p102);

        return dudx.vec(dvdx);
    }

    vec3 tangent_derivative(real pa,real pb,real pc,real u,real v) const {
        // derivative on flat triangle (pa,pb,pc get linearly interpolated):
        // p(u,v) = pa * u + pb * v + pc * (1-u-v)  -> dpdu = pa - pc, dpdv = pb - pc
        // thus dpdx = dpdu*dudx + dpdv*dvdx
        return (pa - pc)*get_dudx(u,v) + (pb - pc)*get_dvdx(u,v);
    }

    vec3 tangent_derivative_at_a(real pa,real pb,real pc) const {
        // note: u = 1, v,w = 0 at a
        vec3 dudx = 3.0*(p300-p201);
        vec3 dvdx = 3.0*(p210-p201);

        vec3 ab(dvdx - dudx);
        vec3 bc(-dvdx);
        
        vec3 n(ab.vec(bc));
        n *= 1.0/n.norm2();

        vec3 grad = (pc - pb)*n.vec(ab)
                  + (pa - pb)*n.vec(bc);
        return grad;
        return (pa - pc)*dudx + (pb - pc)*dvdx;
    }

    vec3 tangent_derivative_at_b(real pa,real pb,real pc) const {
        // note: v = 1, u,w = 0 at b
        vec3 dudx = 3.0*(p120-p021);
        vec3 dvdx = 3.0*(p030-p021);

        vec3 ab(dvdx - dudx);
        vec3 bc(-dvdx);
        
        vec3 n(ab.vec(bc));
        n *= 1.0/n.norm2();

        vec3 grad = (pc - pb)*n.vec(ab)
                  + (pa - pb)*n.vec(bc);
        return grad;
        return (pa - pc)*dudx + (pb - pc)*dvdx;
    }

    vec3 tangent_derivative_at_c(real pa,real pb,real pc) const {
        // note: w = 1, u,v = 0 at c
        vec3 dudx = 3.0*(p102-p003);
        vec3 dvdx = 3.0*(p012-p003);

        vec3 ab(dvdx - dudx);
        vec3 bc(-dvdx);
        
        vec3 n(ab.vec(bc));
        n *= 1.0/n.norm2();

        vec3 grad = (pc - pb)*n.vec(ab)
                  + (pa - pb)*n.vec(bc);
        return grad;
        return (pa - pc)*dudx + (pb - pc)*dvdx;
    }

    vec3 get_normal(real u,real v) const {
        vec3 normal(get_surface_vector(u,v));
        normal.normalize();
        return normal;
    }

    vec3 get_a() const {
        return p300;
    }

    vec3 get_b() const {
        return p030;
    }

    vec3 get_c() const {
        return p003;
    }

private:

    vec3 project_twothrids(vec3 const& p1,vec3 const& p2,vec3 const& n1) const {
        real w = (p1-p2).dot(n1);
        return (2.0*p1 + p2 + w*n1)*(1.0/3.0);
    }

    vec3 p300,p030,p003,p210,p120,p021,p012,p102,p201,p111;
};


} // namespace Bem

#endif // CUBIC_HPP