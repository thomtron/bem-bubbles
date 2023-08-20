#include "ColocSimPin.hpp"
#include "../Integration/Integrator.hpp"
#include "../Mesh/HalfedgeMesh.hpp"
#include "../Mesh/MeshManip.hpp"
#include "../Mesh/MeshIO.hpp"
#include <vector>
#include <numeric> // for iota
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
void ColocSimPin::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
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

    // now we put the terms corresponding to the pinned vertices on the left hand side of the equation G*psi = H*phi
    // The N_pin columns on the right side of G will be overwritten by the negative N_pin columns from H. This is 
    // no problem, since the former would evaluate to zero, since we impose psi=0 on the pinned vertices. The last 
    // N_pin values of x (in G*x = H*phi) will be the potential values at the boundary.

    for(size_t i(m.verts.size()-N_pin);i<m.verts.size();++i) {
        G.col(i) = -H.col(i);
        H.col(i) *= 0.0;
    }

#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif
}


CoordVec ColocSimPin::position_t(Mesh const& m,PotVec& pot) const {
    //cout << "myfunc!" << endl;
    PotVec x;
    CoordVec result = LinLinSim::position_t(m,pot,x);
    for(size_t i(m.verts.size()-N_pin);i<m.verts.size();++i) {
        result[i] = vec3(); // set velocity at pinned nodes to zero
        pot[i] = x[i]; // update the (yet unknown) potential value at the pinned nodes
        //cout << "pot(" << i << ") = " << x[i] << endl;
    }
    return result;
}

void ColocSimPin::rearrange_boundary(Mesh& M,vector<bool> bound) const {
    vector<size_t> permutation(M.verts.size());
    vector<size_t> inverse_permutation(M.verts.size());
    iota(permutation.begin(),permutation.end(),0);
    vector<vec3>& pos(M.verts);
    size_t k(pos.size());
    for(size_t i(0);i<k;++i) {
        if(bound[i]) {
            swap(permutation[i],permutation[k-1]);
            swap(pos[i],pos[k-1]);
            swap(bound[i],bound[k-1]);
            i--;
            k--;
        }
    }
    for(size_t i(0);i<permutation.size();++i) {
        inverse_permutation[permutation[i]] = i;
    }
    for(size_t i(0);i<M.trigs.size();++i) {
        Triplet trig(M.trigs[i]);
        trig.a = inverse_permutation[trig.a];
        trig.b = inverse_permutation[trig.b];
        trig.c = inverse_permutation[trig.c];
        M.trigs[i] = trig;
    }
}

size_t ColocSimPin::set_x_boundary(Mesh& M) const {
    vector<size_t> permutation(M.verts.size());
    vector<size_t> inverse_permutation(M.verts.size());
    iota(permutation.begin(),permutation.end(),0);
    vector<vec3>& pos(M.verts);
    size_t k(pos.size());
    size_t npin(0);
    real threshold = 1e-5;
    for(size_t i(0);i<k;++i) {
        if(abs(pos[i].x)<threshold) {
            pos[i].x = 0.0;
            swap(permutation[i],permutation[k-1]);
            swap(pos[i],pos[k-1]);
            i--;
            k--;
            npin++;
        }
    }
    for(size_t i(0);i<permutation.size();++i) {
        inverse_permutation[permutation[i]] = i;
    }
    for(size_t i(0);i<M.trigs.size();++i) {
        Triplet trig(M.trigs[i]);
        trig.a = inverse_permutation[trig.a];
        trig.b = inverse_permutation[trig.b];
        trig.c = inverse_permutation[trig.c];
        M.trigs[i] = trig;
    }
    cout << "npin = " << npin << endl;
    return npin;
}

void ColocSimPin::remesh(real L) {

    //export_mesh("before_remesh.ply");
    //cout << "           before" << endl;
    
    PotVec new_curv_params = curvature_param();

    curvature_params = damping_factor*curvature_params + (1.0-damping_factor)*new_curv_params;

    if(min_elm_size > 0.0) {
        for(real& elm : curvature_params)
            elm = max(min(elm,1.0/min_elm_size),1.0/max_elm_size);
    }
    
    HalfedgeMesh manip;
    generate_halfedges(manip,mesh);
    
    //if(mesh.check_validity()) cout << "valid." << endl;

    split_edges(manip,curvature_params,L*4.0/3.0);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    //split_edges(manip,curvature_params,L*4.0/3.0);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    //split_edges(manip,curvature_params,L*4.0/3.0);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);

    Mesh new_mesh = generate_mesh(manip);

    // determine the vertices that are on the loop:

    Halfedge* bound = manip.bounds[0];
    Halfedge* u(bound);
    //size_t i(0);
    vector<bool> bounds(manip.verts.size(),false);
    do {
        //cout << u->vert << ", " << i++ << endl;
        bounds[u->vert] = true;
        u = u->next;
    } while(u != bound);

    rearrange_boundary(new_mesh,bounds);
    //cout << (new_mesh.check_validity() ? "valid" : "bad") << endl;

    //if(manip.check_validity()) cout << "halfedgemesh valid. 0" << endl;

    //export_ply("after_collapse_still_valid.ply",new_mesh);

    //cout << "           after_collapse" << endl;

    //size_t num = set_x_boundary(new_mesh); // mixes arrays up!

    assert(bound_inds.size() == N_pin);


    
    CoordVec normals = generate_vertex_normals(new_mesh);

    CoordVec meshverts = new_mesh.verts;
    relax_vertices(new_mesh); // changes vertex positions
    CoordVec pinned(N_pin);
    for(size_t i(0);i<N_pin;++i) {
        pinned[i] = meshverts[meshverts.size()-N_pin+i];
        //cout << pinned[i] << endl;
        new_mesh.verts.pop_back();
        normals.pop_back();
    }

    
    vector<real> new_phi;

    // projecting the new vertices back on the original surface

    project_and_interpolate(new_mesh,normals,new_phi, mesh, make_copy(phi));

    for(vec3 const& elm : pinned) {
        new_mesh.verts.push_back(elm);
        new_phi.push_back(0.0);
    }

    mesh = new_mesh;

    set_phi(new_phi);

    //if(mesh.check_validity()) cout << "valid." << endl;

    /*
    size_t l(0);
    for(vec3 const& elm : mesh.verts) {
        if(abs(elm.x)<1e-5) l++;
    }
    cout << "number of bounds: " << l << endl;

    l = 0;
    for(Triplet const& t : mesh.trigs) {
        size_t nb(0);
        if(abs(mesh.verts[t.a].x)<1e-5) nb++;
        if(abs(mesh.verts[t.b].x)<1e-5) nb++;
        if(abs(mesh.verts[t.c].x)<1e-5) nb++;

        if(nb == 3) throw(out_of_range("daaamn"));
        if(nb == 2) l++;
    }*/

    //cout << "number of bounds (trigs): " << l << endl;

    //export_ply("critical.ply",mesh);
    //cout << "           critical" << endl;

    //generate_halfedges(manip,mesh);
    //if(manip.check_validity()) {
    //    cout << "halfedgemesh valid. 2" << endl;
    //} else {
    //    throw("shit");
    //}
}

/* presumably there is no need to change...
PotVec   ColocSimPin::pot_t(Mesh const& m,CoordVec const& gradients, real t) const {

}

PotVec   ColocSimPin::pot_t_multi(Mesh const& m,CoordVec const& gradients, real t) const {

}*/

} // namespace Bem