#include "LinLinSim.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
#include "../basic/Bem.hpp"
#include "../Mesh/FittingTool.hpp"
#include "../Mesh/Mesh.hpp"
#include "../Mesh/HalfedgeMesh.hpp"
#include "../Mesh/MeshManip.hpp"

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

CoordVec LinLinSim::position_t(Mesh const& m,PotVec& pot) const {
    CoordVec result;

    // setting up the system of equations and solving it.
    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vec3> normals = generate_triangle_normals(m);
    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);
    vector<vec3> tangent_gradients = generate_tangent_gradients(m,pot);

    vector<vec3> vertex_normals = generate_vertex_normals(m);
    vector<vec3> vertex_gradients;

#if LINEAR
    // in the following code block, the gradients at the vertices are determined,
    // by adding the normal derivatives psi_l(i)*normals to the tangential derivatives
    // from generate_tangent_gradients for each triangle and then averaging over the 
    // neighbouring triangles to get the value at the vertex. An alternative would be
    // to locally fit the surface using FittingTool and extract the tangential gradient
    // similarly to the case presented in ConConGalerkinSim.cpp
    for(size_t i(0);i<m.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {

            grad += tangent_gradients[index] + psi_l(i)*normals[index];
            num++;
        }
        grad *= 1.0/num;

        vertex_gradients.push_back(grad);
    }

#else
    // this code does the same thing as the above code, just for bezier triangle interpolation
    // instead of linear interpolation. Here the tangent derivative computation is a bit more complicated, 
    // but we can use a function provided by the Cubic interpolator class.
    vertex_gradients = vector<vec3>(m.verts.size());
    vector<real> num(m.verts.size(),0.0);
    for(size_t i(0);i<m.trigs.size();++i) {
        Triplet t(m.trigs[i]);
        Cubic interp(m.verts[t.a],m.verts[t.b],m.verts[t.c],vertex_normals[t.a],vertex_normals[t.b],vertex_normals[t.c]);
        num[t.a]++; num[t.b]++; num[t.c]++;
        vertex_gradients[t.a] += interp.tangent_derivative_at_a(pot[t.a],pot[t.b],pot[t.c]);
        vertex_gradients[t.b] += interp.tangent_derivative_at_b(pot[t.a],pot[t.b],pot[t.c]);
        vertex_gradients[t.c] += interp.tangent_derivative_at_c(pot[t.a],pot[t.b],pot[t.c]);
    }
    for(size_t i(0);i<m.verts.size();++i) {
        vertex_gradients[i] = vertex_gradients[i]*(1.0/num[i]) + psi_l(i)*vertex_normals[i];
    }
#endif
    result = vertex_gradients;

    return result;

}

PotVec  LinLinSim::pot_t(Mesh const& m,CoordVec const& gradients, real t) const {
    assert(gradients.size() == m.verts.size());
    vector<real> kap(kappa(m));

    real vol(volume(m));

    vector<real> result(gradients.size());
    for(size_t i(0);i<result.size();++i) {
        result[i] = potential_t(gradients[i].norm2(),vol,kap[i],m.verts[i],t);
    }
    return result;
}

PotVec  LinLinSim::pot_t_multi(Mesh const& m, CoordVec const& gradients, real t) const {

    PotVec result(gradients.size());

    vector<vector<size_t>> indices;
    vector<Mesh> meshes = split_by_loose_parts(m,indices);

    vector<vector<vec3>> grads;
    for(vector<size_t> const& inds : indices) {
        vector<vec3> grad;
        for(size_t ind : inds)
            grad.push_back(gradients[ind]);
        grads.push_back(grad);
    }

    for(size_t i(0);i<meshes.size();++i) {
        vector<real> new_pot = pot_t(meshes[i],grads[i],t);
        for(size_t j(0);j<new_pot.size();++j) {
            result[indices[i][j]] = new_pot[j];
        }
    }
    return result;
}

vector<real> LinLinSim::kappa(Mesh const& m) const {
    
    vector<real> kap;
    
    vector<real> gam;
    curvatures(m,kap,gam);
    // 'curvatures' computes the curvature values for each triangle
    // we convert to vertex data by simply averaging

    vector<real> result(m.verts.size(),0.0);
    vector<real> weights(m.verts.size(),0.0);
    for(size_t i(0);i<m.trigs.size();++i) {
        Triplet t(m.trigs[i]);
        result[t.a] += kap[i];
        result[t.b] += kap[i];
        result[t.c] += kap[i];
        weights[t.a]++;
        weights[t.b]++;
        weights[t.c]++;
    }
    for(size_t j(0);j<m.verts.size();++j) {
        result[j]/=weights[j];
    }
    return result;
}

void LinLinSim::evolve_system_RK4(real dp, bool fixdt) {
    // note, position_t is the heavy function here that solves the BEM problem.

    Mesh m2 = mesh;
    Mesh m3 = mesh;
    Mesh m4 = mesh;
    
    CoordVec x1 = mesh.verts;
    PotVec p1 = get_phi();

    CoordVec k1_x = position_t(mesh,p1);
    PotVec k1_p = pot_t_multi(mesh,k1_x,time);

    // check whether using fixed or adaptive timesteps
    real dt;
    if(fixdt) dt = dp;
    else dt = get_dt(dp,k1_x,k1_p);

#ifdef VERBOSE
    cout << "\n dt = " << dt << endl << endl;
#endif
    CoordVec x2 = x1 + (0.5*dt)*k1_x;
    m2.verts = x2;
    PotVec p2 = p1 + (0.5*dt)*k1_p;

    CoordVec k2_x = position_t(m2,p2);
    CoordVec x3 = x1 + (0.5*dt)*k2_x;
    m3.verts = x3;
    PotVec k2_p = pot_t_multi(m2,k2_x,time);
    PotVec p3 = p1 + (0.5*dt)*k2_p;

    CoordVec k3_x = position_t(m3,p3);
    CoordVec x4 = x1 + dt*k3_x;
    m4.verts = x4;
    PotVec k3_p = pot_t_multi(m3,k3_x,time);
    PotVec p4 = p1 + dt*k3_p;

    CoordVec k4_x = position_t(m4,p4);
    
    CoordVec average = (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x)*(1.0/6.0);

    // set a value for psi for export and for guess in solver (this is optional)
    vector<vec3> normals = generate_vertex_normals(mesh);
    psi = Eigen::VectorXd(average.size());
    for(size_t i(0);i<average.size();++i) {
        psi(i) = average[i].dot(normals[i]);
    }

    CoordVec xf = x1 + dt*average;

    PotVec pf = p1 + (dt/6.0)*(pot_t_multi(mesh,k1_x,time) + 2.0*pot_t_multi(m2,k2_x,time) + 2.0*pot_t_multi(m3,k3_x,time) + pot_t_multi(m4,k4_x,time));
    mesh.verts = xf;
    set_phi(pf);
    time += dt;
}

// evolving the system in time with the Euler method. Alternatively evolve_system_RK4 can be used.
void LinLinSim::evolve_system(real dp, bool fixdt) {

    PotVec p = make_copy(phi);

    CoordVec grads = position_t(mesh,p);

    PotVec pot_derivative = pot_t_multi(mesh,grads,time);

    // check whether using fixed or adaptive timesteps
    real dt;
    if(fixdt) dt = dp;
    else dt = get_dt(dp,grads,pot_derivative);
    
#ifdef VERBOSE
    cout << "\n dt = " << dt << endl << endl;
#endif

    mesh.verts = mesh.verts + grads*dt;
    set_phi(p + dt*pot_derivative);
    
    time += dt;

}

// Computes the tangent gradient of phi for each triangle for the case of linear interpolation.
// A vector of vec3's with length m.trigs.size() is returned.
vector<vec3> LinLinSim::generate_tangent_gradients(Mesh const& m, vector<real> const& pot) const {
    vector<vec3> gradients;
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

vector<real> LinLinSim::curvature_param() const {
    vector<real> max_curv = max_curvature(mesh);

    // the following pieces of code can optionally be introduced in order
    // to make the remeshing parameter not only dependent on the maximum 
    // curvature but also on the gradients of phi on the mesh. Since we 
    // haven't found yet the best fitting method, we do not include it in
    // this version. Such adaptions are still subject of experimentation.
    /*
    vector<vec3> tangrad = generate_tangent_gradients(mesh,make_copy(phi));

    vector<vector<size_t>> trig_inds = generate_triangle_indices(mesh);
    vector<real> meandiffgrad(mesh.verts.size());
    for(size_t i(0);i<trig_inds.size();++i){
        real mean = 0.0;
        real num = 0.0;
        for(size_t j : trig_inds[i]){
            for(size_t k : trig_inds[i]){
                if(j != k) {
                    mean += (tangrad[j]-tangrad[k]).norm();
                    num++;
                }
            }
        }
        meandiffgrad[i] = mean/num;
    }


    vector<real> maxgrad(mesh.verts.size(),0.0);
    for(size_t i(0);i<mesh.trigs.size();++i){
        Triplet t(mesh.trigs[i]);
        maxgrad[t.a] = max(maxgrad[t.a],tangrad[i].norm());
        maxgrad[t.b] = max(maxgrad[t.b],tangrad[i].norm());
        maxgrad[t.c] = max(maxgrad[t.c],tangrad[i].norm());
    }

    
    for(size_t i(0);i<max_curv.size();++i) {
        max_curv[i] += 1.0*meandiffgrad[i] + 1.0*maxgrad[i];
    }
    */

    // smoothing the parameters by averaging 
    // over the 2 ring neighbours
    vector<vector<size_t>> two_ring(generate_2_ring(mesh));
    vector<real> max_curv_tmp = max_curv;
    for(size_t i(0);i<mesh.verts.size();++i){
        real mean_curvature = 0.0;
        real num = 0.0;
        for(size_t j : two_ring[i]) {
            mean_curvature += max_curv[j];
            num++;
        }
        max_curv_tmp[i] = mean_curvature/num;
    }
    max_curv = max_curv_tmp;

    return max_curv;
}

void LinLinSim::remesh(real L) {
    
    PotVec new_curv_params = curvature_param();

    curvature_params = damping_factor*curvature_params + (1.0-damping_factor)*new_curv_params;

    if(min_elm_size > 0.0) {
        for(real& elm : curvature_params)
            elm = max(min(elm,1.0/min_elm_size),1.0/max_elm_size);
    }
    
    HalfedgeMesh manip(generate_halfedges(mesh));

    // we iterate here the split_edges, collapse_edges, flip_edges
    // and relax_vertices functions a few times, beginning with one 
    // split edges and then further only using algorithms that keep
    // the triangle count constant or reduce it. Like that the growth 
    // rate of triangle count is intrinsically limited (it cannot
    // explode by just calling remesh once)
    // The application of collapse_edges becomes less costly as 
    // soon as most of the edges verify the condition.

    //split_edges(manip,curvature_params,L*4.0/3.0); 
    
    split_edges(manip,curvature_params,L*0.75); // This "shakes" the mesh up a bit, since most of the
    flip_edges(manip,1);                        // edges are splitted -> gives better stability in some cases
    flip_edges(manip,1);
    relax_vertices(manip);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    relax_vertices(manip);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    relax_vertices(manip);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    relax_vertices(manip);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    relax_vertices(manip);
    flip_edges(manip,1);
    relax_vertices(manip);
    vector<real> new_phi;
    
    Mesh new_mesh = generate_mesh(manip);
    // projecting the new vertices back on the original surface
    project_and_interpolate(new_mesh,new_phi, mesh, make_copy(phi));
    mesh = new_mesh;
    set_phi(new_phi);
}

PotVec compute_exterior_pot(CoordVec const& pos,Mesh M,PotVec phi,PotVec psi) {

    size_t N = pos.size();

    PotVec phi_ext(N);

    #pragma omp parallel
    {
    
    #pragma omp critical 
    {
    cout << "thread " << omp_get_thread_num() << " working..." << endl;
    }

    Integrator inter;
    inter.set_quadrature(quadrature_19);

    #pragma omp for
    for(size_t i = 0;i < N; ++i) {
        Bem::real val(0.0);
        for(Triplet t : M.trigs) {
            val += inter.get_exterior_potential(M.verts,t,phi,psi,pos[i]);
        }
        phi_ext[i] = val;

        if(omp_get_thread_num() == 0)
            cout << " Computing external potential... progress (approx.): " << float(i+1)/N*100.0*omp_get_num_threads() << "%                                    \r" << flush;
    }

    #pragma omp critical 
    {
    cout << "thread " << omp_get_thread_num() << " done." << endl;
    }

    }

    return phi_ext;
}

PotVec LinLinSim::exterior_pot(CoordVec const& positions) const {
    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,mesh);
    Eigen::VectorXd psi_l = solve_system(G,H*phi);

    return compute_exterior_pot(positions,mesh,make_copy(phi),make_copy(psi_l));
}



} // namespace Bem