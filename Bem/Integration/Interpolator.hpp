#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include "../basic/Bem.hpp"

namespace Bem {

// The Interpolator provides functions that return the position of any point of a
// given triangle. It must be initialized with three vectors designating the vertices
// of the triangle in right-handed order.

class Interpolator {
public:

    Interpolator(vec3 a,vec3 b,vec3 c)
        :a(a),ab(b-a),bc(c-b) {
            // corresponds to right-handed normal w.r.t. the order of the given vertices.
            normal_ = ab.vec(bc); 
            area_ = normal_.norm();
            normal_*= (1.0/area_);
        }

    real inline area() const {
        return area_;
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
    real area_;
    
};

} // namespace Bem

#endif // INTERPOLATOR_HPP