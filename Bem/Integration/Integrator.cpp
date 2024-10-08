#include <vector>
#include "Integrator.hpp"
#include "quadrature.hpp"

using Eigen::MatrixXd;

namespace Bem {

// only one specialization here -> this desingularisation is a bit less general b.c. of some analytical results.
// Note further that the identical H factor isn't added here, since it is here a property that is related to the
// local geometry of the mesh and is better added further down.
template<>
void Integrator::integrate_identical_coloc<LinElm>(Interpolator tri_y,LinElm& result) const {
    // Here the case x = tri_y.interpolate(0,0) is handled. 
    // An according reordering of indices happens outside this function.
    // the precision of quad_1d should be chosen higher here than for the galerkin integration

    // M_PI_4 is the variable change factor of transformation from (0,1) to (0,pi/4)
    real jac_factor(0.5*tri_y.area()*M_PI_4);

    LinElm temp;

    // using an adapted rule here since interpolate interpolates not the same unit triangle than in Ning's Work
    
    for(quadrature_1d p : quad_1d){
        p.x *= M_PI_4; // pi/4
        real cosx = cos(p.x);
        real sinx = sin(p.x);
        real dist = (tri_y.interp_relative(cosx,sinx)).norm(); // note: x = a
        real overall_factor = jac_factor*p.weight/dist;
        real factor(1.0/(cosx*cosx));

        // could provide here a function like get_linear_elements()
        temp[0] += overall_factor*cosx*factor;        // a
        temp[1] += overall_factor*(cosx-sinx)*factor; // b
        temp[2] += overall_factor*sinx*factor;        // c
    }

    result = temp;
    return;
}


void Integrator::integrate_identical_coloc_mir(Interpolator tri_y,HomoPair<LinElm>& result) const {
    // Here the case x = tri_y.interpolate(0,0) is handled. 
    // An according reordering of indices happens outside this function.
    // the precision of quad_1d should be chosen higher here than for the galerkin integration

    // M_PI_4 is the variable change factor of transformation from (0,1) to (0,pi/4)
    real jac_factor(0.5*tri_y.area()*M_PI_4);

    LinElm temp;

    // using an adapted rule here since interpolate interpolates not the same unit triangle than in Ning's Work
    
    for(quadrature_1d p : quad_1d){
        p.x *= M_PI_4; // pi/4
        real cosx = cos(p.x);
        real sinx = sin(p.x);
        real dist = (tri_y.interp_relative(cosx,sinx)).norm(); // note: x = a
        real overall_factor = jac_factor*p.weight/dist;
        real factor(1.0/(cosx*cosx));

        // could provide here a function like get_linear_elements()
        temp[0] += overall_factor*cosx*factor;        // a
        temp[1] += overall_factor*(cosx-sinx)*factor; // b
        temp[2] += overall_factor*sinx*factor;        // c
    }

    vec3 a = tri_y.interpolate(0.0,0.0);
    vec3 x = a;
    vec3 b = tri_y.interpolate(1.0,0.0);
    vec3 c = tri_y.interpolate(1.0,1.0);
    a.x = -a.x; // mirror points at x=0 plane
    b.x = -b.x;
    c.x = -c.x;
    Interpolator interp(a,b,c);
    HomoPair<LinElm> disjoint;
    // Important! not mir-version here!
    integrate_disjoint_coloc(x,interp,disjoint); 
    disjoint.H *= (-1.0); // since normal is in 
    // inverse direction due to mirroring!
    disjoint.G += temp;


    result = disjoint;
    
    return;
}


void Integrator::integrate_LinLin(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinLinElm> result;

    integrate(x,tri_i,tri_j,result);
    
    for(size_t i(0);i<3;++i){
        for(size_t j(0);j<3;++j){
            G(tri_i[i],tri_j[j]) += result.G[3*i + j];
            H(tri_i[i],tri_j[j]) += result.H[3*i + j];
        }
    }
    return;
}

void Integrator::integrate_LinLin_local(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    Triplet ref_j = tri_j;

    HomoPair<LinLinElm> result;

    integrate(x,tri_i,tri_j,result);
    
    for(size_t j(0);j<3;++j){

        // handle here the possible permutations on tri_j due to integrate
        size_t ind = 0;
        for(size_t k(0);k<3;++k){
            if(tri_j[j] == ref_j[k]) ind = k;
        }

        for(size_t i(0);i<3;++i){
            G(tri_i[i],ind) += result.G[3*i + j]; // second index on LHS is between 0 and 2 here!
            H(tri_i[i],ind) += result.H[3*i + j];
        }
    }
    
    return;
}

void Integrator::integrate_ConLin(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,MatrixXd& G,MatrixXd& H) const {
    
    Pair<ConElm,LinElm> result;

    integrate(x,tri_i,tri_j,result);

    G(i,j) += result.G;
    for(size_t k(0);k<3;++k) {
        H(i,tri_j[k]) += result.H[k];
    }

    return;
}

void Integrator::integrate_ConLin_local(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,MatrixXd& G,MatrixXd& H) const {
    
    Pair<ConElm,LinElm> result;

    integrate(x,tri_i,tri_j,result);

    G(i,j) += result.G;
    for(size_t k(0);k<3;++k) {
        H(i,k) += result.H[k];
    }

    return;
}

void Integrator::integrate_Con(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,real& G,real& H) const {
    
    HomoPair<real> result;

    integrate(x,tri_i,tri_j,result);
    
    G += result.G;
    H += result.H;
    return;
}

// still the handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
void Integrator::integrate_Lin_coloc_cubic(std::vector<vec3> const& x,std::vector<vec3> const& n,size_t i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinElm> result;

    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {   
        tri_j.cyclic_reorder(i);
        Cubic tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c],n[tri_j.a],n[tri_j.b],n[tri_j.c]);
        integrate_identical_coloc(tri_y,result);
    } else {
        Cubic tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c],n[tri_j.a],n[tri_j.b],n[tri_j.c]);
        integrate_disjoint_coloc(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,tri_j[k]) += result.G[k];
        H(i,tri_j[k]) += result.H[k];
    }
   
}

// The handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
void Integrator::integrate_Lin_coloc(std::vector<vec3> const& x,size_t i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinElm> result;

    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {   
        tri_j.cyclic_reorder(i);
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_identical_coloc(tri_y,result.G); // only G is computed here!
        result.H = 0.0;
    } else {
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_disjoint_coloc(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,tri_j[k]) += result.G[k];
        H(i,tri_j[k]) += result.H[k];
    }
   
}

void Integrator::integrate_Lin_coloc_local_cubic(std::vector<vec3> const& x,std::vector<vec3> const& n,size_t i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinElm> result;

    size_t shift = 0;
    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {  
        if(i == tri_j.b) shift = 1;
        if(i == tri_j.c) shift = 2;
        tri_j.cyclic_reorder(i);
        Cubic tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c],n[tri_j.a],n[tri_j.b],n[tri_j.c]);
        integrate_identical_coloc(tri_y,result);
    } else {
        Cubic tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c],n[tri_j.a],n[tri_j.b],n[tri_j.c]);
        integrate_disjoint_coloc(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,(k+shift)%3) += result.G[k];
        H(i,(k+shift)%3) += result.H[k];
    }
   
}

// The handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
void Integrator::integrate_Lin_coloc_local(std::vector<vec3> const& x,size_t i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinElm> result;

    size_t shift = 0;
    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {  
        if(i == tri_j.b) shift = 1;
        if(i == tri_j.c) shift = 2;
        tri_j.cyclic_reorder(i);
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_identical_coloc(tri_y,result.G); // only G is computed here!
        result.H = 0.0;
    } else {
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_disjoint_coloc(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,(k+shift)%3) += result.G[k];
        H(i,(k+shift)%3) += result.H[k];
    }
   
}

// The handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
void Integrator::integrate_Lin_coloc_local_mir(std::vector<vec3> const& x,size_t i,Triplet tri_j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinElm> result;

    size_t shift = 0;
    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {  
        if(i == tri_j.b) shift = 1;
        if(i == tri_j.c) shift = 2;
        tri_j.cyclic_reorder(i);
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_identical_coloc_mir(tri_y,result); // G and H computed here! - solid angle term still necessary!
    } else {
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        integrate_disjoint_coloc_mir(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,(k+shift)%3) += result.G[k];
        H(i,(k+shift)%3) += result.H[k];
    }
   
}


// Function for computing the potential outside of the mesh surface. x must not be part of the surface!
real Integrator::get_exterior_potential(std::vector<vec3> const& x, Triplet tri_j, std::vector<real> phi, std::vector<real> psi, vec3 y) const {
    HomoPair<LinElm> result;

    Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
    integrate_disjoint_coloc(y,tri_y,result);

    real integral = 
    - result.G[0] * psi[tri_j.a]
    - result.G[1] * psi[tri_j.b]
    - result.G[2] * psi[tri_j.c]

    + result.H[0] * phi[tri_j.a]
    + result.H[1] * phi[tri_j.b]
    + result.H[2] * phi[tri_j.c];

    integral *= (1.0/(4.0*M_PI));  // not necessary for matrix calculation since there the factor cancels!

    return integral;
}



} // namespace Bem
