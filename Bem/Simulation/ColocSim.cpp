#include "ColocSim.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
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

void ColocSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif

    G = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());
    H = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());

    /* // used for completely deterministic runs!
    vector<Eigen::MatrixXd> G_mats(4),H_mats(4);
    for(auto& mat : G_mats)
        mat = G;
    for(auto& mat : H_mats)
        mat = H;
    */

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
            G.col(trip[k]) += G_loc.col(k); //_mats[omp_get_thread_num()]
            H.col(trip[k]) += H_loc.col(k);
        }

        //cout << G << endl << endl << endl;
        }

        
        
    }

    }

/*
    for(auto& mat : G_mats)
        G += mat;
    for(auto& mat : H_mats)
        H += mat;
*/

    // fill the solid angle part of the matrix for cubic interpolation
    for(size_t i(0);i<m.verts.size();++i) {
#if LINEAR
        // fill the solid angle part of the matrix
            
        // 4-pi-rule
        real val_H(0.0);
        //real val_G(0.0);
        for(size_t j(0);j<m.verts.size();++j) {
            val_H -= H(i,j);
            //val_G += G(i,j);
        }
        //cout << val_H << " - " << m.solid_angle_at_vertex(i) << " = " << val_H - m.solid_angle_at_vertex(i) << endl;
        //cout << val_G << " - " << 4.0*M_PI << " = " << val_G - 4.0*M_PI << endl;

        H(i,i) -= (4.0*M_PI - val_H); 
        // note: removed '= -' because for mirror case,
        // H(i,i) already is non zero !

#else
        /*
        real G_row = 0.0;
        real H_row = 0.0;
        for(size_t j(0);j<m.verts.size();++j) {
            G_row += G(i,j);
            H_row += H(i,j);
        }
        cout << G_row << "  " << H_row << "  " << 4*M_PI << "  " << 2*M_PI << endl;
        */
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