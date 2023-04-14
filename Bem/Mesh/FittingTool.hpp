#ifndef FITTINGTOOL_HPP
#define FITTINGTOOL_HPP

#include "../basic/Bem.hpp"

#include <iostream>
#include <vector>
#include <Eigen/Dense>

namespace Bem {

class CoordSystem {
public:

    CoordSystem(vec3 const& origin,vec3 const& x,vec3 y = vec3()) {
        set(origin,x,y);
    }

    CoordSystem() {}


    vec3 transform(vec3 vec) const;
    void set(vec3 const& origin,vec3 const& x,vec3 y = vec3());
    void swap_xz();
    vec3 world_coords(real x, real y, real z) const;
    vec3 world_coords_relative(real x, real y, real z) const;
    vec3 world_coords(vec3 pos) const {
        return world_coords(pos.x,pos.y,pos.z);
    }
    vec3 world_coords_relative(vec3 pos) const {
        return world_coords_relative(pos.x,pos.y,pos.z);
    }

    void print(std::ostream& output) const;

    const vec3& get_X() const {
        return X;
    }
    const vec3& get_Y() const {
        return Y;
    }
    const vec3& get_Z() const {
        return Z;
    }
    const vec3& get_O() const {
        return O;
    }

private:
    vec3 O,X,Y,Z;
};



class FittingTool {
public:
    using real = Bem::real;
    using vec3 = Bem::vec3;

    // fills the vector "fitting_params"
    void compute_quadratic_fit(vec3 normal, vec3 center_vertex, std::vector<vec3> vertices,real d);
    void compute_quadratic_fit(vec3 normal, vec3 center_vertex, std::vector<vec3> vertices);
    vec3 get_position(real x, real y) const;
    real get_curvature() const;
    vec3 get_normal() const;
    std::vector<real> get_params() const {
        return fitting_params;
    }

    CoordSystem copy_coord_system() const {
        return system;
    }

    void print_params() const {
        std::cout << "fitting_params:" << std::endl;
        for(real param : fitting_params) {
            std::cout << param << std::endl;
        }
    }
private:
    std::vector<real> fitting_params;
    CoordSystem system;
};

} // namespace Bem

#endif // FITTINGTOOL_HPP
