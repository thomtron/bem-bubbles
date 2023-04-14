
#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <vector>
#include "../basic/Bem.hpp"
#include "Interpolator.hpp"
#include "Cubic.hpp"
#include "ResultTypes.hpp"

#include "quadrature.hpp"


#include <Eigen/Dense>

namespace Bem {

inline HomoPair<real> integrand(vec3 z,vec3 n);
inline real integrand_identical(vec3 z);

template<typename result_t> inline result_t integrand(real x0,real x1,real y0,real y1,Interpolator interp_x,Interpolator interp_y);
template<typename result_t> inline result_t integrand_identical(real x0,real x1,real y0,real y1,Interpolator interp);

template<typename result_t> inline result_t integrand_coloc(vec3 x,real y0,real y1,Interpolator interp_y);
template<typename result_t> inline result_t integrand_coloc_mir(vec3 x,real y0,real y1,Interpolator interp_y);
template<typename result_t> inline result_t integrand_coloc(vec3 x,real y0,real y1,Cubic const& interp_y);

//template<typename result_t,typename func_type> inline result_t integrand(func_type func,real x0,real x1,Interpolator interp);

class Integrator {
public:

    Integrator()
        :quad_2d(quadrature_3)
        ,quad_1d(gauss_3) {}


    void integrate_LinLin               (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_LinLin_local         (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_disLinLin            (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_ConLin               (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_ConLin_local         (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Con                  (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,real& G,real& H) const;
    void integrate_Lin_coloc_cubic      (std::vector<vec3> const& x,std::vector<vec3> const& n,size_t i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Lin_coloc            (std::vector<vec3> const& x,size_t i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Lin_coloc_local_cubic(std::vector<vec3> const& x,std::vector<vec3> const& n,size_t i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Lin_coloc_local      (std::vector<vec3> const& x,size_t i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Lin_coloc_local_mir  (std::vector<vec3> const& x,size_t i,Triplet tri_j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    void integrate_Lin_discoloc         (std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,Eigen::MatrixXd& G,Eigen::MatrixXd& H) const;
    
    real get_exterior_potential(std::vector<vec3> const& x, Triplet tri_j, std::vector<real> phi, std::vector<real> psi, vec3 y) const;

    // for the following templates, the passed function must be analytic everywhere - 
    // no treatment of singularities is applied.
    template<typename func_type>
    void integrate_function(Interpolator tri,func_type func,real& result) const;
    // func_t must be able to be called with (real,real,Interpolator,element_t) as argument, returning a real
    template<typename func_t,typename element_t>
    void integrate_function(Interpolator tri, func_t func,element_t approx_func,real& result) const;
    // "manual" version of above
    template<typename result_t> 
    void integrate_function(Interpolator tri, result_t (*func)(real,real,Interpolator),result_t& result) const;

    void set_quadrature(std::vector<quadrature_2d> const& quad) {
        quad_2d = quad;
    }

    void set_quadrature(std::vector<quadrature_1d> const& quad) {
        quad_1d = quad;
    }

//private:

    template<typename result_t>
    void integrate_disjoint     (Interpolator tri_x,Interpolator tri_y, result_t& result) const;
    template<typename result_t>
    void integrate_shared_vertex(Interpolator tri_x,Interpolator tri_y, result_t& result) const;
    template<typename result_t>
    void integrate_shared_edge  (Interpolator tri_x,Interpolator tri_y, result_t& result) const;
    template<typename result_t>
    void integrate_identical    (Interpolator tri_x, result_t& result) const;
    //template<typename result_t>
    //void integrate_identical_spec (Interpolator tri_x, result_t& result) const;


    template<typename result_t>
    void integrate_disjoint_coloc (vec3 x,Cubic const& tri_y,result_t& result) const;

    template<typename result_t>
    void integrate_disjoint_coloc (vec3 x,Interpolator tri_y,result_t& result) const;
    
    template<typename result_t>
    void integrate_disjoint_coloc_mir (vec3 x,Interpolator tri_y,result_t& result) const;
    

    template<typename result_t>
    void integrate_identical_coloc (Cubic const& tri_y,result_t& result) const;

    template<typename result_t>
    void integrate_identical_coloc (Interpolator tri_y,result_t& result) const;

    void integrate_identical_coloc_mir (Interpolator tri_y,HomoPair<LinElm>& result) const;


    template<typename result_t>
    void integrate(std::vector<vec3> const& x,Triplet& tri_i,Triplet& tri_j,result_t& result) const;


    std::vector<quadrature_2d> quad_2d;
    std::vector<quadrature_1d> quad_1d;
    
};


// inline functions:

// basic functions
inline HomoPair<real> integrand(vec3 z,vec3 n) {
    //real inv_dist = z.norm2(); // test A -> res: 4/3 resp. 1/3
    //real inv_dist = sin(z.x)*sin(z.y); // test B -> res: 0 resp. 0
    real inv_dist = 1.0/z.norm();
    return HomoPair<real>(inv_dist,-z.dot(n)*inv_dist*inv_dist*inv_dist);
}

inline real integrand_identical(vec3 z) {
    //return z.norm2(); // test A -> res: 4/3 resp. 1/3
    //return sin(z.x)*sin(z.y); // test B -> res: 0 resp. 0
    return 1.0/z.norm();
}

template<>
inline HomoPair<real> integrand<HomoPair<real>>(real x0,real x1,real y0,real y1,Interpolator interp_x,Interpolator interp_y) {
    return integrand(interp_y.interpolate(y0,y1)-interp_x.interpolate(x0,x1),interp_y.normal());
}

template<>
inline HomoPair<LinLinElm> integrand<HomoPair<LinLinElm>>(real x0,real x1,real y0,real y1,Interpolator interp_x,Interpolator interp_y) { 
    HomoPair<real> basic(integrand<HomoPair<real>>(x0,x1,y0,y1,interp_x,interp_y));
    HomoPair<LinLinElm> result;

    result.G = get_linear_elements(x0,x1,y0,y1);
    result.H = result.G;

    result.G *= basic.G;
    result.H *= basic.H;
    return result;
}

template<>
inline real integrand_identical<real>(real x0,real x1,real y0,real y1,Interpolator interp) {
    return integrand_identical(interp.interpolate(y0,y1) - interp.interpolate(x0,x1));
}

template<>
inline LinLinElm integrand_identical<LinLinElm>(real x0,real x1,real y0,real y1,Interpolator interp) {
    LinLinElm result(get_linear_elements(x0,x1,y0,y1));
    result *= integrand_identical<real>(x0,x1,y0,y1,interp);
    return result;
}

template<>
inline Pair<ConElm,LinElm> integrand<Pair<ConElm,LinElm>>(real x0,real x1,real y0,real y1,Interpolator interp_x,Interpolator interp_y) { 
    HomoPair<real> basic(integrand<HomoPair<real>>(x0,x1,y0,y1,interp_x,interp_y));
    Pair<ConElm,LinElm> result;

    result.G = basic.G;
    result.H = get_linear_elements(y0,y1);
    result.H *= basic.H;
    return result;
}


template<>
inline HomoPair<real> integrand_coloc<HomoPair<real>>(vec3 x,real y0,real y1,Cubic const& interp_y) {
    vec3 normal(interp_y.get_surface_vector(y0,y1));
    HomoPair<real> basic(integrand(interp_y.interpolate(y0,y1) - x,normal));
    basic *= normal.norm();
    return basic;
}

template<>
inline HomoPair<real> integrand_coloc<HomoPair<real>>(vec3 x,real y0,real y1,Interpolator interp_y) {
    return integrand(interp_y.interpolate(y0,y1) - x,interp_y.normal());
}

template<>
inline HomoPair<LinElm> integrand_coloc<HomoPair<LinElm>>(vec3 x,real y0,real y1,Cubic const& interp_y) {
    vec3 normal(interp_y.get_surface_vector(y0,y1));
    real jac = normal.norm();
    normal *= (1.0/jac); // normalize here manually, since norm() is already computed
    HomoPair<real> basic(integrand(interp_y.interpolate(y0,y1) - x,normal));
    basic *= jac;
    HomoPair<LinElm> result;

    result.G = get_linear_elements_for_cubic(y0,y1);
    result.H = result.G;
    
    result.G *= basic.G;
    result.H *= basic.H;
    return result;
}

template<>
inline HomoPair<LinElm> integrand_coloc<HomoPair<LinElm>>(vec3 x,real y0,real y1,Interpolator interp_y) {
    HomoPair<real> basic(integrand(interp_y.interpolate(y0,y1) - x,interp_y.normal()));
    
    HomoPair<LinElm> result;

    result.G = get_linear_elements(y0,y1);
    result.H = result.G;
    
    result.G *= basic.G;
    result.H *= basic.H;
    return result;
}

template<>
inline HomoPair<LinElm> integrand_coloc_mir<HomoPair<LinElm>>(vec3 x,real y0,real y1,Interpolator interp_y) {
    // HomoPair<real> basic(integrand(interp_y.interpolate(y0,y1) - x,interp_y.normal()));
    vec3 y = interp_y.interpolate(y0,y1);
    vec3 n = interp_y.normal();
    HomoPair<real> basic(integrand(y - x,n));
    y.x = -y.x; // mirror
    n.x = -n.x;
    basic += integrand(y - x,n);

    HomoPair<LinElm> result;

    result.G = get_linear_elements(y0,y1);
    result.H = result.G;
    
    result.G *= basic.G;
    result.H *= basic.H;

    return result;
}


/*
template<>
inline real integrand<real,real (*)(vec3)>(real (*func)(vec3),real x0,real x1,Interpolator interp) {
    return func(interp.interpolate(x0,x1));
}

template<>
inline real integrand<real,real (*)(vec3,vec3)>(real (*func)(vec3 x,vec3 n),real x0,real x1,Interpolator interp) {
    return func(interp.interpolate(x0,x1),interp.normal());
}

template<>
inline LinElm integrand<LinElm,real (*)(vec3)>(real (*func)(vec3),real x0,real x1,Interpolator interp) {
    LinElm result(get_linear_elements(x0,x1));
    result *= func(interp.interpolate(x0,x1));
    return result;
}

template<>
inline LinElm integrand<LinElm,real (*)(vec3 x,vec3 n)>(real (*func)(vec3 x,vec3 n),real x0,real x1,Interpolator interp) {
    LinElm result(get_linear_elements(x0,x1));
    result *= func(interp.interpolate(x0,x1),interp.normal());
    return result;
}
*/

// implementations of member templates:


// galerkin methods
template<typename result_t>
void Integrator::integrate_disjoint(Interpolator tri_x,Interpolator tri_y, result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    for(quadrature_2d const& q_i : quad_2d) {
        for(quadrature_2d const& q_j : quad_2d) {
            result_t temp(integrand<result_t>(
                                q_i.x+q_i.y,q_i.y, // the addition of q_i/j.y is necessary to transform the quadrature to 
                                q_j.x+q_j.y,q_j.y, // the other unit triangle with 45° corners at (0,0) and (1,1)
                                tri_x,
                                tri_y
                                ));
                                
            temp *= q_i.weight*q_j.weight;

            result += temp;
        }
    }

    result *= tri_x.jacobian()*tri_y.jacobian();

    return;
}

template<typename result_t>
void Integrator::integrate_identical(Interpolator tri_x, result_t& result) const {
    result.G = 0.0;
    result.H = result_t::identical_H_factor;
    result.H *= -tri_x.jacobian(); // negative factor here!

    using G_t = typename result_t::G_t;

    for(quadrature_1d xi : quad_1d) {
        for(quadrature_1d eta1 : quad_1d) {
            for(quadrature_1d eta2 : quad_1d) {
                for(quadrature_1d eta3 : quad_1d) {

                    real weight = xi.weight
                               *eta1.weight
                               *eta2.weight
                               *eta3.weight
                               // and jacobian of duffy coords:
                               *xi.x*xi.x*xi.x 
                               *eta1.x*eta1.x
                               *eta2.x;


                    real A = xi.x;
                    real B = A*eta1.x;
                    real C = B*eta2.x;
                    real D = C*eta3.x;

                    // result.H is zero here: normal is orthogonal to y-x

                    G_t temp_G(0.0);
                    
                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,A-B+C
                                    ,A-D
                                    ,A-B
                                    ,tri_x);
                               

                    temp_G += integrand_identical<G_t>(  
                                     A-D
                                    ,A-B
                                    ,A
                                    ,A-B+C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,B-C+D
                                    ,A-C
                                    ,B-C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A-C
                                    ,B-C
                                    ,A
                                    ,B-C+D
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A-D
                                    ,B-D
                                    ,A
                                    ,B-C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,B-C
                                    ,A-D
                                    ,B-D
                                    ,tri_x);

                    temp_G *= weight;

                    result.G += temp_G;
                }
            }
        }
    }

    result.G *= tri_x.jacobian()*tri_x.jacobian();

    return;
}

/*
template<>
void Integrator::integrate_identical_spec(Interpolator tri_x, HomoPair<LinElm>& result) const {
    result.G = 0.0;
    result.H = result_t::identical_H_factor;
    result.H *= tri_x.jacobian();

    using G_t = typename result_t::G_t;

    for(quadrature_1d xi : quad_1d) {
        for(quadrature_1d eta1 : quad_1d) {
            for(quadrature_1d eta2 : quad_1d) {
                for(quadrature_1d eta3 : quad_1d) {

                    real weight = xi.weight
                               *eta1.weight
                               *eta2.weight
                               *eta3.weight
                               // and jacobian of duffy coords:
                               *xi.x*xi.x*xi.x 
                               *eta1.x*eta1.x
                               *eta2.x;


                    real A = xi.x;
                    real B = A*eta1.x;
                    real C = B*eta2.x;
                    real D = C*eta3.x;

                    // result.H is zero here: normal is orthogonal to y-x

                    G_t temp_G(0.0);
                    
                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,A-B+C
                                    ,A-D
                                    ,A-B
                                    ,tri_x);
                               

                    temp_G += integrand_identical<G_t>(  
                                     A-D
                                    ,A-B
                                    ,A
                                    ,A-B+C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,B-C+D
                                    ,A-C
                                    ,B-C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A-C
                                    ,B-C
                                    ,A
                                    ,B-C+D
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A-D
                                    ,B-D
                                    ,A
                                    ,B-C
                                    ,tri_x);

                    temp_G += integrand_identical<G_t>(  
                                     A
                                    ,B-C
                                    ,A-D
                                    ,B-D
                                    ,tri_x);

                    temp_G *= weight;

                    result.G += temp_G;
                }
            }
        }
    }

    result.G *= tri_x.jacobian()*tri_x.jacobian();

    return;
}
*/

template<typename result_t>
void Integrator::integrate_shared_edge(Interpolator tri_x,Interpolator tri_y, result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    for(quadrature_1d xi : quad_1d) {
        for(quadrature_1d eta1 : quad_1d) {
            for(quadrature_1d eta2 : quad_1d) {
                for(quadrature_1d eta3 : quad_1d) {

                    real weight = xi.weight
                               *eta1.weight
                               *eta2.weight
                               *eta3.weight
                               // and jacobian of duffy coords:
                               *xi.x*xi.x*xi.x
                               *eta1.x*eta1.x;
                               // attention! below we add eta2.x to the jacobian for some of the terms

                    real A = xi.x;        // xi*...
                    real B = A*eta1.x;    // eta1
                    real C = B*eta2.x;    // eta1*eta2
                    real D = C*eta3.x;    // eta1*eta2*eta3

                    result_t temp;

                    temp += integrand<result_t>(
                                  A
                                 ,B
                                 ,A-D
                                 ,C-D
                                 ,tri_x
                                 ,tri_y);

                    temp += integrand<result_t>(
                                  A-C
                                 ,B-C
                                 ,A
                                 ,D
                                 ,tri_x
                                 ,tri_y);
                    
                    temp += integrand<result_t>(
                                  A-D
                                 ,C-D
                                 ,A
                                 ,B
                                 ,tri_x
                                 ,tri_y);
                    
                    temp += integrand<result_t>(
                                  A-D
                                 ,B-D
                                 ,A
                                 ,C
                                 ,tri_x
                                 ,tri_y);
                    
                    temp *= eta2.x;

                    temp += integrand<result_t>(
                                  A
                                 ,B*eta3.x // special case here
                                 ,A-C
                                 ,B-C
                                 ,tri_x
                                 ,tri_y);

                    temp *= weight;

                    result += temp;

                }
            }
        }
    }

    result *= tri_x.jacobian()*tri_y.jacobian();

    return;
}

template<typename result_t>
void Integrator::integrate_shared_vertex(Interpolator tri_x,Interpolator tri_y, result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    for(quadrature_1d xi : quad_1d) {
        for(quadrature_1d eta1 : quad_1d) {
            for(quadrature_1d eta2 : quad_1d) {
                for(quadrature_1d eta3 : quad_1d) {

                    real weight = xi.weight
                               *eta1.weight
                               *eta2.weight
                               *eta3.weight
                               // and jacobian of duffy coords:
                               *xi.x*xi.x*xi.x
                               *eta2.x;

                    real A = xi.x;
                    real B = A*eta1.x;
                    real C = A*eta2.x; // different to the other cases!
                    real D = C*eta3.x;

                    result_t temp;
                    temp += integrand<result_t>(A,B,C,D,tri_x,tri_y);
                    temp += integrand<result_t>(C,D,A,B,tri_x,tri_y);
                    temp *= weight;

                    result += temp;

                }
            }
        }
    }

    result *= tri_x.jacobian()*tri_y.jacobian();

    return;
}

// colocation methods
template<typename result_t>
void Integrator::integrate_disjoint_coloc(vec3 x,Interpolator tri_y,result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    //QuadratureList_2d quad = *triangle_quadratures[3];

    for(quadrature_2d q_y : quad_2d) {
        result_t temp(integrand_coloc<result_t>(
                            x,
                            q_y.x+q_y.y,q_y.y, // transform to other unit triangle
                            tri_y));
                            
        temp *= q_y.weight;

        result += temp;
    }
    

    result *= tri_y.jacobian();

    return;
}

template<typename result_t>
void Integrator::integrate_disjoint_coloc_mir(vec3 x,Interpolator tri_y,result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    //QuadratureList_2d quad = *triangle_quadratures[3];

    for(quadrature_2d q_y : quad_2d) {
        result_t temp(integrand_coloc_mir<result_t>(
                            x,
                            q_y.x+q_y.y,q_y.y, // transform to other unit triangle
                            tri_y));
                            
        temp *= q_y.weight;

        result += temp;
    }
    

    result *= tri_y.jacobian();
    
    return;
}

template<typename result_t>
void Integrator::integrate_disjoint_coloc(vec3 x,Cubic const& tri_y,result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;

    //QuadratureList_2d quad = *triangle_quadratures[3];

    for(quadrature_2d q_y : quad_2d) {
        result_t temp(integrand_coloc<result_t>(
                            x,
                            q_y.x,q_y.y,
                            tri_y));
                            
        temp *= q_y.weight;

        result += temp;
    }

    return;
}

template<typename result_t>
void Integrator::integrate_identical_coloc(Cubic const& tri_y,result_t& result) const {
    result.G = 0.0;
    result.H = 0.0;
    vec3 x = tri_y.get_a();



    // I use here an easy version of duffy coords: e1 = u, e2 = u*v, u,v in (0,1) (gauss. quadr.)
    // The integration factor of this transform is u.
    
    for(quadrature_1d p : quad_1d){

        for(quadrature_1d q : quad_1d){

            real u = p.x;
            real v = q.x;
            real e1 = u;
            real e2 = u*v;

            result_t temp(integrand_coloc<result_t>(
                            x,
                            1.0-e1,e2, // flip triangle such that it fits the triangle representation of Cubic
                            tri_y));
                            
            temp *= p.weight*q.weight*u;

            result += temp;
        }
    }
    
}



// function integration

// expects a function func that accepts a position vector and a normal vector as arguments.
template<typename func_type>
void Integrator::integrate_function(Interpolator tri,func_type func,real& result) const {
    result = 0.0;

    for(quadrature_2d const& q_y : quad_2d) {
        real temp(func(tri.interpolate(q_y.x+q_y.y,q_y.y) // transform to other unit triangle
                          ,tri.normal()));
                            
        temp *= q_y.weight;

        result += temp;
    }

    result *= tri.jacobian();

    return;
}

template<typename result_t>
void Integrator::integrate_function(Interpolator tri, result_t (*func)(real x0,real x1,Interpolator interp),result_t& result) const {
    result = 0.0;

    for(quadrature_2d const& q_y : quad_2d) {
        result_t temp(func( q_y.x+q_y.y,q_y.y, // transform to other unit triangle
                            tri));
                            
        temp *= q_y.weight;

        result += temp;
    }

    result *= tri.jacobian();

    return;
}

template<typename func_t,typename element_t>
void Integrator::integrate_function(Interpolator tri, func_t func,element_t approx_func,real& result) const {
    result = 0.0;

    for(quadrature_2d const& q_y : quad_2d) {
        real temp(func( q_y.x+q_y.y,q_y.y, // transform to other unit triangle
                        tri,
                        approx_func));
                            
        temp *= q_y.weight;

        result += temp;
    }

    result *= tri.jacobian();

    return;
}

/*
template<typename func_type>
void Integrator::get_linear_representation(func_type func,std::vector<vec3> const& x,Triplet tri,std::vector<real>& result) const {
    Interpolator interp(x[tri.a],x[tri.b],x[tri.c]);
    LinElm res;
    integrate_function(interp,func,res);
    result[tri.a] = res[0];
    result[tri.b] = res[1];
    result[tri.c] = res[2];
    return;
}
*/




template<typename result_t>
void Integrator::integrate(std::vector<vec3> const& x,Triplet& tri_i,Triplet& tri_j,result_t& result) const {

    size_t num_equals(0);             
    for(size_t i :     {tri_i.a,tri_i.b,tri_i.c}){
        for(size_t j : {tri_j.a,tri_j.b,tri_j.c}){
            if(i == j) {
                num_equals++;
            }
        }
    }

    switch(num_equals){

        case 0: // disjoint triangles
        {
            Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);
            Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);

            integrate_disjoint(tri_x,tri_y,result);
            //result.G = 0.0;
            break;
        }
        case 1: // shared vertex
        {

            for(size_t j(0);j<3;++j){
                for(size_t i(0);i<3;++i){
                    if(tri_i[i] == tri_j[j]) {
                        tri_i.cyclic_reorder(tri_i[i]);
                        tri_j.cyclic_reorder(tri_j[j]);
                        // exit the loops:
                        i = 3;
                        j = 3;
                    }
                }
            }

            Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);
            Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);

            integrate_shared_vertex(tri_x,tri_y,result);
            //result.G = 0.0;
            break;
        }
        case 2: // shared edge
        {
            for(size_t j(0);j<3;++j){
                for(size_t i(0);i<3;++i){
                    if(tri_i[i] == tri_j[j]) {
                        tri_i.cyclic_reorder(tri_i[i]);
                        tri_j.cyclic_reorder(tri_j[j]);
                        // exit the loops:
                        i = 3;
                        j = 3;
                    }
                }
            }

            real sign(1.0);
            
            std::swap(tri_i.b,tri_i.c);
            if(tri_i.b != tri_j.b) {
                std::swap(tri_i.b,tri_i.c);
                std::swap(tri_j.b,tri_j.c);
                sign = -1.0;  // -> minus sign due to normal pointing in other direction 
                              //    in the y-integration (corresponding to t_j). Note that 
                              //    for the triangle t_i and for the G integrals it makes 
                              //    no difference, since no dot product with the normal has
                              //    to be evaluated! 
            }

            Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);
            Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);

            //result_t testbit;
            integrate_shared_edge(tri_x,tri_y,result);
            result.H *= sign;
            //testbit.G = result.G;
            //integrate_shared_edge(tri_y,tri_x,result);
            //result.G -= testbit.G;
            break;
        }
        case 3: // identical triangle
        {
            Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);

            integrate_identical(tri_x,result);
            //result.G = 0.0;
            break;
        }
        default:
            throw "damn, there is an invalid index!";
    }

    return;
}

} // namespace Bem

#endif // INTEGRATOR_HPP
