#ifndef RESULTTYPES_HPP
#define RESULTTYPES_HPP

#include <array>
#include "../basic/Bem.hpp"


namespace Bem {

// note for more complicated statements it would make sense to 
// create lazy evaluation methods for ElementArray, but I think
// this isn't so useful in our case.
template<size_t N>
class ElementArray : public std::array<real,N>{
public:
    ElementArray(real const& value) {
        std::array<real,N>::fill(value);
    }
    ElementArray() 
        :ElementArray(0.0) {}
    ElementArray(std::array<real,N> array)
        :std::array<real,N>(array) {}
    
    void operator*=(real const& value) {
        for(real& elm : *this){
            elm *= value;
        }
    }
    void operator+=(ElementArray const& other) {
        for(size_t i(0);i<N;++i){
            (*this)[i] += other[i];
        }
    }
    void operator-=(ElementArray const& other) {
        for(size_t i(0);i<N;++i){
            (*this)[i] -= other[i];
        }
    }
    void operator=(real const& value) {
        std::array<real,N>::fill(value);
    }
};

template<typename A, typename B>
struct Pair {
    using G_t = A;
    using H_t = B;
    static const H_t identical_H_factor;

    G_t G;
    H_t H;
    Pair(real g,real h)
        :G(g),H(h) {}

    Pair()
        :Pair(0.0,0.0) {}

    void operator+=(Pair const& other) {
        G += other.G;
        H += other.H;
    }

    void operator-=(Pair const& other) {
        G -= other.G;
        H -= other.H;
    }

    void operator*=(real const& scalar) {
        G *= scalar;
        H *= scalar;
    }
};

template<typename T>
using HomoPair = Pair<T,T>;

using LinElm = ElementArray<3>;
using LinLinElm = ElementArray<9>;
using ConElm = real;

inline LinLinElm get_linear_elements(real x0,real x1,real y0,real y1);
inline LinElm get_linear_elements(real y0,real y1);
// same for ConstantLinear 

inline LinLinElm get_linear_elements(real x0,real x1,real y0,real y1) {
    real x0y0(x0*y0),x0y1(x0*y1),x1y0(x1*y0),x1y1(x1*y1);

    LinLinElm functions;
    // triangles with coordinates:          x - y
    functions[0] = (1.0-x0-y0+x0y0);     // a - a
    functions[1] = (y0-y1-x0y0+x0y1);    // a - b
    functions[2] = (y1-x0y1);            // a - c
    functions[3] = (x0-x1-x0y0+x1y0);    // b - a
    functions[4] = (x0y0-x0y1-x1y0+x1y1);// b - b
    functions[5] = (x0y1-x1y1);          // b - c
    functions[6] = (x1-x1y0);            // c - a
    functions[7] = (x1y0-x1y1);          // c - b
    functions[8] = (x1y1);               // c - c

    return functions;
}

inline LinElm get_linear_elements(real x0,real x1) {
    LinElm functions;
    functions[0] = 1.0 - x0; // a
    functions[1] = x0 - x1;  // b
    functions[2] = x1;       // c
    return functions;
}

inline LinElm get_linear_elements_for_cubic(real x0,real x1) {
    LinElm functions;
    functions[0] = x0;            // a
    functions[1] = x1;            // b
    functions[2] = 1.0 - x0 - x1; // c
    return functions;
}

} // namespace Bem


#endif // RESULTTYPES_HPP