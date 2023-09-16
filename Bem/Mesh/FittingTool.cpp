#include "FittingTool.hpp"
#include <cassert>
#include <iostream>

using namespace std;

namespace Bem {

void CoordSystem::set(vec3 const& origin,vec3 const& x,vec3 y) {
    O = origin;
    X = x;
    assert(!X.null(1e-10) && "x is zero!"); // x must not be zero

    X.normalize();

    if(y.null()) {
        if(X.y == 0.0) {
            y = vec3(0.0,1.0,0.0);
        } else {
            y = vec3(0.0,-X.z/X.y,1.0);
        }
    }

    assert(!X.vec(y).null(1e-7) && "x and y are collinear!"); // vector x and y must not be collinear

    Y = y - X.dot(y)*X;
    Y.normalize();

    Z = X.vec(Y);
}

vec3 CoordSystem::transform(vec3 vec) const {
    vec -= O;
    return vec3(vec.dot(X),vec.dot(Y),vec.dot(Z));
}

// permute axis cyclically such that the given X corresponds now to the Z axis
void CoordSystem::swap_xz() {
    swap(X,Z);
    swap(X,Y);
}

void CoordSystem::print(ostream& output) const {
    output << "O: " << O << endl;
    output << "X: " << X << endl;
    output << "Y: " << Y << endl;
    output << "Z: " << Z << endl;
}

vec3 CoordSystem::world_coords(real x, real y, real z) const {
    return vec3(O + x*X + y*Y + z*Z);
}

vec3 CoordSystem::world_coords_relative(real x, real y, real z) const {
    return vec3(x*X + y*Y + z*Z);
}

void FittingTool::compute_quadratic_fit(vec3 normal, vec3 center_vertex, std::vector<vec3> vertices) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    system.set(center_vertex,normal);
    system.swap_xz();

    real d(0.0);
    for(vec3& vec : vertices) {
        vec = system.transform(vec);
        d += vec.norm();
    }
    d /= double(vertices.size());

    size_t n(vertices.size());
    if(n<6) {
        //if(n<5) cerr << "Danger: less than 5 points were given to FittingTool!" << endl;
        //else cerr << "Warning: less than 6 points were given to FittingTool!" << endl;
    }
    MatrixXd B(n,6);
    MatrixXd W = MatrixXd::Zero(n,n);
    VectorXd zeta(n);
    for(size_t i(0);i<n;++i){
        vec3 const& vec(vertices[i]);
        
        B(i,0) = 1.0;
        B(i,1) = vec.x;
        B(i,2) = vec.y;
        B(i,3) = vec.x*vec.y;
        B(i,4) = vec.x*vec.x;
        B(i,5) = vec.y*vec.y;

        zeta(i) = vec.z;
        // the weight of each point depends exponentially on its distance to the origin
        W(i,i) = exp(-vec.norm()/(2*d));
    }

    MatrixXd A = B.transpose()*W*B;
    VectorXd b = B.transpose()*W*zeta;

    Eigen::FullPivLU<Eigen::MatrixXd> solver;
    solver.compute(A);
    VectorXd x = solver.solve(b);
    
    fitting_params.clear();
    for(size_t i(0);i<6;++i){
        fitting_params.push_back(x(i));
    }
}

// returns the woorld coordinates of a point on the fit surface at (x,y)
vec3 FittingTool::get_position(real x, real y) const {
    return system.world_coords(x,y,fitting_params[0]
                                  +fitting_params[1]*x
                                  +fitting_params[2]*y
                                  +fitting_params[3]*x*y
                                  +fitting_params[4]*x*x
                                  +fitting_params[5]*y*y);
}

// see Wang_2014
Bem::real FittingTool::get_curvature() const {
    real a1,a2,a3,a4,a5;
    a1 = fitting_params[1];
    a2 = fitting_params[2];
    a3 = fitting_params[3];
    a4 = fitting_params[4];
    a5 = fitting_params[5];

    return -(a4 + a5 + a5*a1*a1 + a4*a2*a2 - a1*a2*a3)/pow((1.0 + a1*a1 + a2*a2),1.5);
}

vec3 FittingTool::get_normal() const {
    // compute gradient of phi(x,y,z) := f_a(x,y)-z at (x,y,z) = (0,0,0)
    vec3 grad(fitting_params[1],fitting_params[2],-1.0);
    grad.normalize();
    return -system.world_coords_relative(grad.x,grad.y,grad.z);
}

} // namespace Bem