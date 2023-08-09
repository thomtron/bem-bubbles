#ifndef BEM_HPP
#define BEM_HPP

#include "vector3d.hpp"
#include "Triplet.hpp"
#include <vector>
#include <cassert>

namespace Bem {
    using real = double;
    using vec3 = vector3d<real>;

    // two typedefs that are useful for representing scalar functions 
    // and vector fields on the bubble surface
    using CoordVec = std::vector<vec3>;
    using PotVec = std::vector<real>;

    // The following templates enable basic arithmetic with std-vectors of values
    template<typename T>
    std::vector<T> operator*(real s,std::vector<T> vec) {
        for(T& elm : vec)
            elm *= s;
        return vec;
    }
    template<typename T>
    std::vector<T> operator*(std::vector<T> vec, real s) {
        return s*vec;
    }
    template<typename T>
    std::vector<T> operator+(std::vector<T> const& v1,std::vector<T> v2) {
        assert(v1.size() == v2.size());
        for(size_t i(0);i<v1.size();++i) {
            v2[i] += v1[i];
        }
        return v2;
    }
}

#endif // BEM_HPP
