#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include "../basic/Bem.hpp"

namespace Bem {

class Interpolator {
public:

    Interpolator(vec3 a,vec3 b,vec3 c)
        :a(a),ab(b-a),bc(c-b) {
            // corresponds to right-handed normal w.r.t. order of given vertices.
            normal_ = ab.vec(bc); 
            jacobian_ = normal_.norm();
            normal_*= (1.0/jacobian_);
        }

    real inline jacobian() const {
        return jacobian_;
    }

    vec3 inline normal() const {
        return normal_;
    }

    vec3 inline interpolate(real u,real v) const {
        return a + u*ab + v*bc;
    }

    vec3 inline interp_relative(real u,real v) const {
        return u*ab + v*bc;
    }

private:
    vec3 a,ab,bc,normal_;
    real jacobian_;
    
};

} // namespace Bem

#endif // INTERPOLATOR_HPP