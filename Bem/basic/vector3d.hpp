#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <cmath> // for sqrt()

// vector3d is a lightweight 3D vector class, that supports basic vector arithmetic such as
// addition, multiplication by scalars and dot- and vector product.Further one can normalize
// vectors and test whether they are zero (optionally with a tolerance parameter).

// the type real must support all operations that a floating point number supports.
template<typename real>
class vector3d {
public:
    vector3d(real x,real y,real z)
        :x(x),y(y),z(z) {}
    // default initialization with 0,0,0
    vector3d() :vector3d(0,0,0) {}
    // (default copy and destructor)

    // I make them public to avoid making getters and setters
    real x;
    real y;
    real z;

    vector3d& operator+=(vector3d const& other) {
        x+=other.x;
        y+=other.y;
        z+=other.z;
        return *this;
    }

    vector3d& operator-=(vector3d const& other) {
        x-=other.x;
        y-=other.y;
        z-=other.z;
        return *this;
    }

    vector3d& operator*=(real scalar) {
        x*=scalar;
        y*=scalar;
        z*=scalar;
        return *this;
    }

    vector3d operator-() const;

    real dot(vector3d const& other) const {
        return x*other.x + y*other.y + z*other.z;
    }

    vector3d vec(vector3d const& other) const;

    real norm2() const {
        return dot(*this);
    }

    real norm() const {
        return sqrt(norm2());
    }

    void normalize() {
        *this *= 1.0/norm();
    }

    bool null(real epsilon = 0.0) {
        return std::abs(x) <= epsilon and std::abs(y) <= epsilon and std::abs(z) <= epsilon;
    }

};

// member functions:

template<typename real> 
vector3d<real> vector3d<real>::vec(vector3d<real> const& other) const {
    vector3d result;
    result.x = y*other.z - z*other.y;
    result.y = z*other.x - x*other.z;
    result.z = x*other.y - y*other.x;

    return result;
}

template<typename real>
vector3d<real> vector3d<real>::operator-() const {
    return *this*(-1.0);
}

// non member functions:

template<typename real> 
vector3d<real> operator*(vector3d<real> v,real s) {
    v *= s;
    return v;
}

template<typename real> 
vector3d<real> operator*(real s,vector3d<real> v) {
    v *= s;
    return v;
}
template<typename real> 
vector3d<real> operator+(vector3d<real> v1,vector3d<real> const& v2) {
    v1 += v2;
    return v1;
}
template<typename real> 
vector3d<real> operator-(vector3d<real> v1,vector3d<real> const& v2) {
    v1 -= v2;
    return v1;
}

#include <iostream>
#include <iomanip>
template<typename real> 
std::ostream& operator<<(std::ostream& output,vector3d<real> const& v) {
#ifdef VERBOSE
    if( (abs(v.x) > 1e-4 or v.x == 0.0) and (abs(v.y) > 1e-4 or v.y == 0.0) and (abs(v.z) > 1e-4 or v.z == 0.0) ) {
        output << std::fixed << std::setprecision(4);
    } else {
        output << std::scientific << std::setprecision(4);
    }
    output << "( " << v.x << " , " << v.y << " , " << v.z << " )" << std::defaultfloat;
#else
    output << v.x << ' ' << v.y << ' ' << v.z << ' ';
#endif
    return output;
}

#endif // VECTOR3D_HPP
