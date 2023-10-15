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
        Eigen::VectorXd tmp(G.col(i));
        G.col(i) = -H.col(i);
        H.col(i) = -tmp;
    }

#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif
}

const real beta = 0.1;
//constexpr real costheta_eq = cos(106.0/180.0*M_PI); // = cos(106°)
const real costheta_eq = cos(130.0/180.0*M_PI); // = cos(106°)


// costheta = normal * (1,0,0) = normal.x (angle taken from interior of buble to interface
// and normal is surf normal at contact line)
real drift_velocity(real costheta) {
    return beta*(costheta_eq - costheta);
}


CoordVec ColocSimPin::position_t(Mesh const& m,PotVec& pot) const {

    CoordVec result;

    // pot contains phi values for tangent_gradients
    vector<vec3> tangent_gradients = generate_tangent_gradients(m,pot);

    // but the N_pin last elements must be filled with the normal gradients
    // before solving the system (note that for RK4, they may be updated in the function...)
    for(size_t i(m.verts.size()-N_pin);i<m.verts.size();++i) {
        pot[i] = psi(i);
    }

    // setting up the system of equations and solving it.
    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vec3> normals = generate_triangle_normals(m);
    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);

    vector<vec3> vertex_normals = generate_vertex_normals(m);
    vector<vec3> vertex_gradients;

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

    result = vertex_gradients;

    for(size_t i(m.verts.size()-N_pin);i<m.verts.size();++i) {
        vec3 normal(vertex_normals[i]);
        real costheta = normal.x;
        normal.x = 0.0;
        normal.normalize(); // this vector is tangential to the solid surface  
        result[i] = normal*drift_velocity(costheta);

        //result[i] = vec3(); // set velocity at pinned nodes to zero
        

        pot[i] = psi_l(i); // update the (yet unknown) potential value at the pinned nodes
        
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
    
    HalfedgeMesh manip(mesh);

    // copy doesn't work properly right now!! there are still some errors to be solved!
    HalfedgeMesh orig_h(mesh);
    vector<size_t> perm;
    Halfedge* start = orig_h.bounds[0];
    Halfedge* curr = start;
    do {
        vec3 pos = orig_h.vpos[curr->vert];
        size_t ind_min(mesh.verts.size()-N_pin);
        real dist_min((mesh.verts[ind_min]-pos).norm2());
        for(size_t j(1);j<N_pin;++j) {
            size_t l(mesh.verts.size()-N_pin+j);
            real dist = (mesh.verts[l]-pos).norm2();
            if(dist < dist_min) {
                dist_min = dist;
                ind_min = l;
            }
        }

        perm.push_back(ind_min);
        
        curr = curr->next;

    } while (curr != start);
    //cout << "did i survive so far?a" << endl;
    
    //if(mesh.check_validity()) cout << "valid." << endl;

    //export_ply("before_split.ply",generate_mesh(manip));
    
    split_edges(manip,curvature_params,L*4.0/3.0);
    //export_ply("after_split.ply",generate_mesh(manip));
    //if(manip.check_validity()) cout << "valid." << endl;
    //generate_halfedges(manip,generate_mesh(manip));
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    //export_ply("last_valid.ply",generate_mesh(manip));
    //if(manip.check_validity()) cout << "valid." << endl;
    //else throw;
    flip_edges(manip,1);
    flip_edges(manip,1);
    //split_edges(manip,curvature_params,L*4.0/3.0);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    //export_ply("prob_invalid.ply",generate_mesh(manip));
    //if(manip.check_validity()) cout << "valid." << endl;
    flip_edges(manip,1);
    flip_edges(manip,1);
    //if(manip.check_validity()) cout << "valid." << endl;
    //split_edges(manip,curvature_params,L*4.0/3.0);
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);

    //cout << "did i survive so far?ab" << endl;

    Mesh new_mesh = generate_mesh(manip);

    //cout << "did i survive so far?b" << endl;

    // determine the vertices that are on the loop:

    Halfedge* bound = manip.bounds[0];
    Halfedge* u(bound);
    size_t n_b(0);
    vector<bool> bounds(manip.verts.size(),false);
    do {
        //cout << u->vert << ", " << i++ << endl;
        bounds[u->vert] = true;
        u = u->next;
        n_b++;
    } while(u != bound);

    rearrange_boundary(new_mesh,bounds);

    Mesh new_mesh_II(new_mesh);
    //cout << (new_mesh.check_validity() ? "valid" : "bad") << endl;

    //if(manip.check_validity()) cout << "halfedgemesh valid. 0" << endl;

    //export_ply("after_collapse_still_valid.ply",new_mesh);

    //cout << "           after_collapse" << endl;

    //size_t num = set_x_boundary(new_mesh); // mixes arrays up!

    N_pin = n_b;

    //cout << "did i survive so far?c" << endl;
    
    CoordVec normals = generate_vertex_normals(new_mesh);

    CoordVec pinned,pinnednormals;
    for(size_t i(0);i<N_pin;++i) {
        pinned.push_back(new_mesh.verts[new_mesh.verts.size()-N_pin+i]);
        vec3 normal = normals[new_mesh.verts.size()-N_pin+i];
        normal.x = 0.0;
        normal.normalize();
        pinnednormals.push_back(normal);
    }

    relax_vertices(new_mesh); // changes vertex positions

    // IMPORTANT NOTE: we assume here that the boundary of the mesh remains unchanged!
    // dammit they get permuted of course !
    for(size_t i(0);i<N_pin;++i) {
        new_mesh.verts.pop_back();
        normals.pop_back();
    }

    
    vector<real> new_phi,new_psi;

    // projecting the new vertices back on the original surface

    project_and_interpolate(new_mesh,normals,new_phi,new_psi, mesh, make_copy(phi), make_copy(psi));

    //cout << "did i survive so far?d" << endl;


    vector<vec3> manip_norms = generate_vertex_normals(manip);
    vector<vec3> better_norms(N_pin);

    Halfedge* ha = bound;

    vec3 old2(manip.vpos[ha->vert]);
    ha = ha->next;
    vec3 old(manip.vpos[ha->vert]);
    ha = ha->next;
    vec3 cur(manip.vpos[ha->vert]);

    Halfedge* hb = ha;
    do {

        // get permutation:

        vec3 pos = manip.vpos[ha->vert];
        size_t ind_min(0);
        real dist_min((new_mesh_II.verts[new_mesh_II.verts.size()-N_pin]-pos).norm2());
        for(size_t j(1);j<N_pin;++j) {
            size_t l(new_mesh_II.verts.size()-N_pin+j);
            real dist = (new_mesh_II.verts[l]-pos).norm2();
            if(dist < dist_min) {
                dist_min = dist;
                ind_min = j;
            }
        }
        //cout << ind_min << "  " << dist_min << endl;

        // done

        vec3 nor = manip_norms[ha->vert];
        vec3 nor_orig = nor;
        nor.x = 0.0;
        nor.normalize();

        ha = ha->next;
        old2 = old;
        old = cur;
        cur = manip.vpos[ha->vert];

        vec3 nor_a = (cur - old).vec(vec3(1,0,0));
        vec3 nor_b = (old - old2).vec(vec3(1,0,0));
        nor_a.normalize();
        nor_b.normalize();

        vec3 nor_o = nor_a + nor_b;
        nor_o.normalize();

        if(nor_o.dot(nor_orig)<0.0) {
            nor_o = -nor_o;
        }

        better_norms[ind_min] = nor_o;

        real val = (nor_o - nor).norm();
        if(val > 0.1) {
            cout << "val = " << val << ", nor = " << nor << ", nor_o = " << nor_o << endl; 
        }
        
    } while(ha != hb);

    Mesh borderprojects;
    size_t trigind_offset = 0;
    PotVec borderpots;

    for(size_t i(0);i<N_pin;++i) {

        // project the pinned ones too on to the old ring
        vec3 pos(pinned[i]);
        //vec3 nor(pinnednormals[i]);
        vec3 nor(better_norms[i]);

        /*

        vec3 nor_a = (new_mesh.verts[perm[(i+1)%kp]] - new_mesh.verts[perm[i]]).vec(vec3(1,0,0));
        vec3 nor_b = (new_mesh.verts[perm[i]] - new_mesh.verts[perm[(i-1)%kp]]).vec(vec3(1,0,0));
        nor_a.normalize();
        nor_b.normalize();

        vec3 nor_o = nor_a + nor_b;
        nor_o.normalize();

        cout << (nor_o - nor).norm() << ", nor = " << nor << ", nor_o = " << nor_o << endl;
        */

        real s_min = 0.0;
        real t_min = 0.0;
        size_t ind_min = 0;
        
        bool first(true);
        size_t kp(perm.size());
        for(size_t j(0);j<kp;j++) {
            vec3 k = mesh.verts[perm[(j+1)%kp]] - mesh.verts[perm[j]]; // boundary element vector
            vec3 O = mesh.verts[perm[j]] - pos;

            real t = (O.vec(nor)).x/(nor.vec(k)).x;
            real s = nor.dot(O + t*k)/nor.norm2();

            if(t >= 0.0 and t < 1.0) {
                // ray is intersecting current line segment of boundary
                if(first or (abs(s) < abs(s_min))) {
                    first = false;
                    s_min = s;
                    t_min = t;
                    ind_min = j;
                }
            }
        }
        if(first) throw(out_of_range("project_and_interpolate: no intersection at boundary"));


        vec3 pos_a = mesh.verts[perm[ind_min]];
        vec3 pos_b = mesh.verts[perm[(ind_min+1)%kp]];

        vec3 old_pos = pos - nor*0.2;

        pos += s_min*nor;

        borderprojects.verts.push_back(old_pos);
        borderprojects.verts.push_back(pos_a);
        borderprojects.verts.push_back(pos);
        borderprojects.verts.push_back(pos_b);

        borderprojects.trigs.push_back(Triplet(trigind_offset,trigind_offset+1,trigind_offset+2));
        borderprojects.trigs.push_back(Triplet(trigind_offset,trigind_offset+2,trigind_offset+3));
        trigind_offset += 4;

        //t_min = 1.0 - t_min;

        real phi_int = t_min*phi(perm[ind_min]) + (1.0-t_min)*phi(perm[(ind_min+1)%kp]);
        real psi_int = t_min*psi(perm[ind_min]) + (1.0-t_min)*psi(perm[(ind_min+1)%kp]);

        borderpots.push_back(0.0);
        borderpots.push_back(phi(perm[ind_min]));
        borderpots.push_back(phi_int);
        borderpots.push_back(phi(perm[(ind_min+1)%kp]));

        
        

        new_mesh.verts.push_back(pos);
        new_phi.push_back(phi_int);
        new_psi.push_back(psi_int);
    }

    borderprojects.add(mesh); // add the old mesh!
    for(auto const& elm : make_copy(phi)) {
        borderpots.push_back(elm);
    }

    export_ply_float("borderprojects/mesh-"+to_string(index++)+".ply",borderprojects,borderpots);

    mesh = new_mesh;

    set_phi(new_phi);
    set_psi(new_psi);

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
