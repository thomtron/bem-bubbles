#include "GalerkinSim.hpp"
#include "../Integration/Integrator.hpp"
#include <vector>
#include <omp.h>


#ifdef VERBOSE
#include <chrono>
#endif

using namespace std;
#ifdef VERBOSE
using namespace chrono;
#endif

namespace Bem {

// assemble_matrices computes the matrix elements of the system matrices G and H for the given Mesh m.
void GalerkinSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif

    G = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());
    H = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
    
    Mesh local(m);
    const vector<vec3>& x(local.verts);
    const vector<vec3> n(generate_vertex_normals(local));
    size_t N(local.verts.size());
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
        Eigen::MatrixXd G_loc = Eigen::MatrixXd::Zero(N,3);
        Eigen::MatrixXd H_loc = Eigen::MatrixXd::Zero(N,3);
        const Triplet trip(local.trigs[j]);
        
        for(size_t i(0);i<M;++i) {
            int_local.integrate_LinLin_local(x,local.trigs[i],trip,G_loc,H_loc);
        }

#ifdef VERBOSE
        if(omp_get_thread_num() == 0)
            cout << " Assembling matrices... progress (approx.): " << float(j+1)/M*100.0*omp_get_num_threads() << "%                                    \r" << flush;
#endif
        #pragma omp critical 
        {
        for(size_t k(0);k<3;++k) {
            G.col(trip[k]) += G_loc.col(k);
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
}

} // namespace Bem