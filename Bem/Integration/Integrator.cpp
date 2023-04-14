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

    // M_PI_4 is the jacobian of transformation from (0,1) to (0,pi/4)
    real jac_factor(0.5*tri_y.jacobian()*M_PI_4);

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

    // M_PI_4 is the jacobian of transformation from (0,1) to (0,pi/4)
    real jac_factor(0.5*tri_y.jacobian()*M_PI_4);

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

#include <iostream>
using namespace std;

void Integrator::integrate_disLinLin(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,MatrixXd& G,MatrixXd& H) const {
    
    HomoPair<LinLinElm> result;

    // handle equal stuff! (tri_i/tri_j transformations)

    Triplet tri_i_ref(tri_i);
    Triplet tri_j_ref(tri_j);

    integrate(x,tri_i,tri_j,result);

    std::vector<size_t> ind_i(3),ind_j(3);
    for(size_t k(0);k<3;++k){
        for(size_t l(0);l<3;++l){
            if(tri_i[k] == tri_i_ref[l]){
                ind_i[k] = l;
                l = 3; // exit loop
            } 
        }
    }
    for(size_t k(0);k<3;++k){
        for(size_t l(0);l<3;++l){
            if(tri_j[k] == tri_j_ref[l]){
                ind_j[k] = l;
                l = 3; // exit loop
            } 
        }
    }

    //if(tri_i_ref[ind_i[0]] != tri_i[0]) cout << "fail 0." << endl;
    //else cout << "nofail." << endl;

    //if(tri_i[ind_i[0]] != tri_i_ref[0]) cout << ind_i[0] << '.' << ind_i[1] << '.' << ind_i[2] << '-' << tri_i[0] << '.' << tri_i[1] << '.' << tri_i[2] << '-' << tri_i_ref[0] << '.' << tri_i_ref[1] << '.' << tri_i_ref[2] << endl;
    
    for(size_t k(0);k<3;++k){
        for(size_t l(0);l<3;++l){
            G(3*i + ind_i[k],3*j + ind_j[l])  += result.G[3*k + l];
            H(3*i + ind_i[k],tri_j[l])        += result.H[3*k + l];
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
        //Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        Cubic tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c],n[tri_j.a],n[tri_j.b],n[tri_j.c]);
        integrate_disjoint_coloc(x[i],tri_y,result);
    }

    for(size_t k(0);k<3;++k) {
        G(i,tri_j[k]) += result.G[k];
        H(i,tri_j[k]) += result.H[k];
    }
   
}

// still the handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
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

// still the handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information

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

// still the handling of the solid angle factor is left for a function in Simulation which has access to all the mesh information
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

void Integrator::integrate_Lin_discoloc(std::vector<vec3> const& x,Triplet tri_i,Triplet tri_j,size_t i,size_t j,MatrixXd& G,MatrixXd& H) const {
    
    std::vector<HomoPair<LinElm>> result(3);

    if(i == j) {
        // handle same triangle
        Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);
        
        for(size_t k(0);k<3;++k){
            quadrature_2d q(quadrature_3[k]);

            LinElm temp;
            vec3 x_coloc = tri_x.interpolate(q.x+q.y,q.y);

            Interpolator   tri_y(x_coloc,x[tri_j.a],x[tri_j.b]);
            integrate_identical_coloc(tri_y,temp);
            result[k].G += temp;

            tri_y = Interpolator(x_coloc,x[tri_j.b],x[tri_j.c]);
            integrate_identical_coloc(tri_y,temp);
            result[k].G += temp;

            tri_y = Interpolator(x_coloc,x[tri_j.c],x[tri_j.a]);
            integrate_identical_coloc(tri_y,temp);
            result[k].G += temp;
            
        }
    } else {
        // handle disjoint triangles
        Interpolator tri_x(x[tri_i.a],x[tri_i.b],x[tri_i.c]);
        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
        
        for(size_t k(0);k<3;++k){
            quadrature_2d q(quadrature_3[k]);
            vec3 x_coloc = tri_x.interpolate(q.x+q.y,q.y);

            integrate_disjoint_coloc(x_coloc,tri_y,result[k]);
        }

    }


//    if(i == tri_j.a or i == tri_j.b or i == tri_j.c) {   
//        if(i == tri_j.b) {
//            indices.cyclic_reorder(1);
//        } else if(i == tri_j.c) {
//            indices.cyclic_reorder(2);
//        }
//        tri_j.cyclic_reorder(i);
//        
//        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
//        integrate_identical_coloc(tri_y,result.G); // only G is computed here!
//        result.H = 0.0;
//    } else {
//        Interpolator tri_y(x[tri_j.a],x[tri_j.b],x[tri_j.c]);
//        integrate_disjoint_coloc(x[i],tri_y,result);
//    }


    
    for(size_t l(0);l<3;++l){
        for(size_t k(0);k<3;++k){
            G(3*i + l,3*j + k)  += result[l].G[k];
            H(3*i + l,tri_j[k]) += result[l].H[k];
        }
    }

    // Also: entweder 3 pönkt pro drüüegg auz kolokationspönkt nää (am beschte nemmsch quadraturpönkt - dänk ech - 
    // het zwor de nooteil, dass d'nodepönkt säuber ned denne send...) oder du mosches met de lineare elemänt mache.
    // lineari galerkin elemänt chiemed gloub ohni zuesätzlechi rächnige bem matrize uufstöue uus

    // wemmer done 3 kolokationspönkt wetti nää, mösstme d'integration vom singulare drüüegg aapasse! -> chame mache,
    // e dem mer s'drüegg i 3 chlineri omwandlet, wo ei egge jewius am kolokationsponkt het. das chönntme tatsächlech 
    // no ennerhaub vo dere fonktion done mache. denn aber statt en index, es drüüegg för i entgägenää!
   
}

// x must not be part of the boundary surface
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

    //std::cout << result.G[0] << ' ' << phi[tri_j.a] << endl;
    return integral;
}



} // namespace Bem
