#include "ColocSim.hpp"
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

// assemble_matrices computes the matrix elements of the system matrices G and H for the given 
// Mesh m. We discern two times two cases: linear interpolation or cubic interpolation (LINEAR) and
// Mirrored mesh or not (for the "image charges" technique) (MIRROR_MESH)
void ColocSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif

    G = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());
    H = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
    
    Mesh local(m);
    const CoordVec& x(local.verts);
    const CoordVec n(generate_vertex_normals(local));
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
        
        for(size_t i(0);i<N;++i) {
#if LINEAR
    #if MIRROR_MESH
            int_local.integrate_Lin_coloc_local_mir(x,i,trip,G_loc,H_loc);
    #else
            int_local.integrate_Lin_coloc_local(x,i,trip,G_loc,H_loc);

    #endif
#else
            int_local.integrate_Lin_coloc_local_cubic(x,n,i,trip,G_loc,H_loc);
#endif
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

    // The following part fills the solid angle depending part of the matrix.
    // We use here the summation of the H-matrix elements. alternatively, the solid angle
    // can be evaluated with solid_angle_at_vertex(Mesh,size_t) to prevent the summation
    for(size_t i(0);i<m.verts.size();++i) {
#if LINEAR
            
        // '4-pi-rule'
        real val_H(0.0);
        
        for(size_t j(0);j<m.verts.size();++j) {
            val_H -= H(i,j);
        }
        //alternative: val_H = solid_angle_at_vertex(m,i);

        H(i,i) -= (4.0*M_PI - val_H); 

#else
        // if cubic bezier triangle interpolation is used, the solid angle always is 1/2 at the vertices
        H(i,i) -= 2.0*M_PI;
#endif
    }


#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif
}

} // namespace Bem