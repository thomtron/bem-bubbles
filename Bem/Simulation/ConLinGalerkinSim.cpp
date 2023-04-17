#include "ConLinGalerkinSim.hpp"
#include "../Integration/Integrator.hpp"
#include "../basic/Bem.hpp"
#include <vector>
#ifdef VERBOSE
#include <chrono>
#endif

using namespace std;
#ifdef VERBOSE
using namespace chrono;
#endif
using namespace Bem;



void ConLinGalerkinSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif
    G = Eigen::MatrixXd::Zero(m.trigs.size(),m.trigs.size());
    H = Eigen::MatrixXd::Zero(m.trigs.size(),m.verts.size());


    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
    
    Mesh local(m);
    const vector<vec3>& x(local.verts);
    size_t M(local.trigs.size());
    Integrator int_local(inter);

#ifdef VERBOSE
    #pragma omp master
    {
        cout << "number of threads: " << omp_get_num_threads() << endl;
        cout << "number of cpu's:   " << omp_get_num_procs() << endl;
    }
#endif
    
    #pragma omp for
    for(size_t j = 0;j<M;++j) {
        Eigen::MatrixXd H_loc = Eigen::MatrixXd::Zero(M,3);
        const Triplet trip(local.trigs[j]);
        
        for(size_t i(0);i<M;++i) {
            int_local.integrate_ConLin_local(x,local.trigs[i],trip,i,j,G,H_loc);
        }

#ifdef VERBOSE
        if(omp_get_thread_num() == 0)
            cout << " Assembling matrices... progress (approx.): " << float(j+1)/M*100.0*omp_get_num_threads() << "%                                    \r" << flush;
#endif
        #pragma omp critical 
        {
        for(size_t k(0);k<3;++k) {
            H.col(trip[k]) += H_loc.col(k);
        }
        }
  
    }
    }

#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif

    return;
}

CoordVec ConLinGalerkinSim::position_t(Mesh const& m,PotVec const& pot) const {

    // setting up the system of equations and solving it.
    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);
    vector<vec3> normals = generate_triangle_normals(m);
    vector<vec3> tangent_gradients = generate_tangent_gradients(m,pot);
    vector<vec3> vertex_gradients;

    // adding the tangent derivatives of phi to the normal derivatives as in LinLinSim.cpp
    for(size_t i(0);i<m.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {

            grad += tangent_gradients[index] + psi_l(index)*normals[index];
            num++;
        }
        grad *= 1.0/num;

        vertex_gradients.push_back(grad);
    }

    return vertex_gradients;

}

// same function as in LinLinSim.cpp
CoordVec ConLinGalerkinSim::generate_tangent_gradients(Mesh const& m, PotVec const& pot) const {
    CoordVec gradients;
    for(Triplet trig : m.trigs) {
        // this formula for the gradient of linear functions on 
        // an unstructured grid is taken from: Dombre_2019, page 8.
        vec3 ab(m.verts[trig.b]-m.verts[trig.a]);
        vec3 bc(m.verts[trig.c]-m.verts[trig.b]);
        
        vec3 n(ab.vec(bc));
        n *= 1.0/n.norm2();

        vec3 grad = (pot[trig.c]-pot[trig.b])*n.vec(ab)
                  + (pot[trig.a]-pot[trig.b])*n.vec(bc);

        gradients.push_back(grad);
    }
    return gradients;
}
