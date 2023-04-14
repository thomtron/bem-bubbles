#include "ColocSim-group.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
#include "../basic/Bem.hpp"
#include "../Mesh/FittingTool.hpp"
#include "../Mesh/HalfedgeMesh.hpp"
#include "../Mesh/MeshManip.hpp"
#include "../Mesh/MeshIO.hpp"
//#include "../Mesh/CurvatureTensor.hpp"
#include "SimData.hpp"
#include <vector>
#include <omp.h>

/*
    Notes:
    Most stable results were obtained until now with fitting the polynomial only
    the direct neighbours of the vertex. large sigma leads to instability quite 
    soon. but with a small sigma too the simulation becomes unstable after one or
    two oscillations...

    I found out that the curvatures obtained with the formula are always larger than 1...

*/

#define ANALYTICS false

#define LINEAR false

#ifdef VERBOSE
#include <chrono>
#endif

using namespace std;
#ifdef VERBOSE
using namespace chrono;
#endif

namespace Bem {

namespace Group {

void ColocSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif

    G = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());
    H = Eigen::MatrixXd::Zero(m.verts.size(),m.verts.size());

    omp_set_num_threads(100);

    #pragma omp parallel
    {

    
    Mesh local(m);
    const vector<vec3>& x(local.verts);
    const vector<vec3> n(generate_vertex_normals(local));
    size_t N(local.verts.size());
    size_t M(local.trigs.size());

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
            inter.integrate_Lin_coloc_local(x,i,trip,G_loc,H_loc);
#else
            inter.integrate_Lin_coloc_local_cubic(x,n,i,trip,G_loc,H_loc);
#endif
        }

#ifdef VERBOSE
        if(omp_get_thread_num() == 0)
            cout << " Assembling matrices... progress (approx.): " << float(j+1)/M*100.0*omp_get_num_threads() << "%                                    \r" << flush;
#endif
        #pragma omp critical 
        for(size_t k(0);k<3;++k) {
            G.col(trip[k]) += G_loc.col(k);
            H.col(trip[k]) += H_loc.col(k);
        }

        //cout << G << endl << endl;
        
    }

    }

    // fill the solid angle part of the matrix for cubic interpolation
    for(size_t i(0);i<m.verts.size();++i) {
#if LINEAR
        // fill the solid angle part of the matrix
            
        // 4-pi-rule
        real val_H(0.0);
        //real val_G(0.0);
        for(size_t j(0);j<m.verts.size();++j) {
            val_H += H(i,j);
            //val_G += G(i,j);
        }
        //cout << val_H << " - " << m.solid_angle_at_vertex(i) << " = " << val_H - m.solid_angle_at_vertex(i) << endl;
        //cout << val_G << " - " << 4.0*M_PI << " = " << val_G - 4.0*M_PI << endl;

        H(i,i) = + (4.0*M_PI - val_H);

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
        H(i,i) += 2.0*M_PI;
#endif
    }

    

#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif
}

/*
SimData ColocSim::gradient(SimData const& X) const {
    // make a new mesh object, with the positions given in X
    // copy the triangles vector from "mesh" into that mesh.

    Mesh m;
    m.verts = X.pos;
    m.trigs = mesh.trigs;

    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi = solve_system(G,H*make_copy(X.phi));

    SimData result;

    vector<real> kappa;
    vector<vec3> vertex_gradients = generate_gradients(kappa,m,make_copy(X.phi),psi);


    for(size_t i(0);i<m.verts.size();++i) {

        result.pos.push_back(vertex_gradients[i]);
        result.phi.push_back(potential_t(vertex_gradients[i].norm2(),kappa[i]));
    }

    return result;
   
    / *

    // test here where phi is the radial velocity

    Mesh m;
    m.verts = X.pos;
    m.trigs = mesh.trigs;

    SimData result;
    result.pos = m.generate_vertex_normals_max();
    result.phi = vector<real>(result.pos.size());

    real R(0.0);
    for(vec3 const& x : X.pos)
        R += x.norm();
    R /= static_cast<real>(X.pos.size());

    for(size_t i(0);i<result.pos.size();++i) {
        result.pos[i] *= X.phi[i];

        result.phi[i] = -1.0/R*(3.0/2.0* X.phi[i]*X.phi[i] + 2.0*sigma/R + 1.0 - epsilon*pow(1.0/R,3*lambda));
    }

    return result;

    * /
}
*/


CoordVec ColocSim::position_t(Mesh const& m,PotVec const& pot) const {
    CoordVec result;

    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vec3> normals = generate_triangle_normals(m);
    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);
    vector<vec3> tangent_gradients = generate_tangent_gradients(m,pot);
    /*
    for(size_t i(0);i<m.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {
            grad += tangent_gradients[index]+psi_l(i)*normals[index];
            num++;
        }
        grad *= 1.0/num;

        result.push_back(grad);
    }*/

    vector<vec3> vertex_normals = generate_vertex_normals(m);
    vector<vec3> vertex_gradients;


#if LINEAR
    for(size_t i(0);i<m.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {

            grad += tangent_gradients[index] + psi_l(i)*normals[index];
            //grad += triangle_gradients[index];
            num++;
        }
        grad *= 1.0/num;

        //grad -= grad.dot(vertex_normals[i])*vertex_normals[i]; // remove normal part
        //grad += vertex_normals[i]*psi(i);

        vertex_gradients.push_back(grad);

        //cout << grad.norm() << endl;
    }
#else
    vertex_gradients = vector<vec3>(m.verts.size());
    vector<real> num(m.verts.size(),0.0);
    for(size_t i(0);i<m.trigs.size();++i) {
        Triplet t(m.trigs[i]);
        Cubic interp(m.verts[t.a],m.verts[t.b],m.verts[t.c],vertex_normals[t.a],vertex_normals[t.b],vertex_normals[t.c]);
        num[t.a]++; num[t.b]++; num[t.c]++;
        vertex_gradients[t.a] += interp.tangent_derivative_at_a(pot[t.a],pot[t.b],pot[t.c]);
        vertex_gradients[t.b] += interp.tangent_derivative_at_b(pot[t.a],pot[t.b],pot[t.c]);
        vertex_gradients[t.c] += interp.tangent_derivative_at_c(pot[t.a],pot[t.b],pot[t.c]);
/*
        cout << interp.tangent_derivative_at_a(1.0,0.0,0.0) << endl;
        cout << interp.tangent_derivative_at_b(1.0,0.0,0.0) << endl;
        cout << interp.tangent_derivative_at_c(1.0,0.0,0.0) << endl;
        cout << interp.tangent_derivative_at_a(pot[t.a],pot[t.b],pot[t.c]) << endl;
        cout << interp.tangent_derivative_at_b(pot[t.a],pot[t.b],pot[t.c]) << endl;
        cout << interp.tangent_derivative_at_c(pot[t.a],pot[t.b],pot[t.c]) << endl;
        cout << tangent_gradients[i] << endl << endl;
        */
    }
    for(size_t i(0);i<m.verts.size();++i) {
        vertex_gradients[i] = vertex_gradients[i]*(1.0/num[i]) + psi_l(i)*vertex_normals[i];
    }
#endif

    result = vertex_gradients;

    return result;

}

PotVec  ColocSim::pot_t(Mesh const& m,CoordVec const& gradients) const {
    assert(gradients.size() == m.verts.size());
    vector<real> kap(kappa(m));
    real vol(volume(m));

    vector<real> result(gradients.size());
    for(size_t i(0);i<result.size();++i) {
        result[i] = potential_t(gradients[i].norm2(),vol,kap[i]);
    }
    return result;
}

vector<real> ColocSim::kappa(Mesh const& m) const {
    
    vector<real> kap;
    /*
    vector<vector<size_t>> verts = generate_neighbours(m);

    for(size_t i(0);i<verts.size();++i) {
        verts[i].push_back(i);
    }
    
    
    vector<vec3> norms = generate_vertex_normals(m);

    size_t n(m.verts.size());

    for(size_t i(0);i<n;++i) {
        vector<vec3> positions;
        for(size_t j : verts[i]) {
            positions.push_back(m.verts[j]);
        }
        FittingTool fit;
        fit.compute_quadratic_fit(norms[i],m.verts[i],positions);
        kap.push_back(fit.get_curvature());
    }*/
    /*
    Bem::real mean_rad(0.0);
    for(vec3 x : m.verts) {
        mean_rad += x.norm();
    }
    mean_rad /= static_cast<Bem::real>(m.verts.size());
    for(size_t i(0);i<m.verts.size();++i) {
        kap.push_back(1.0/mean_rad);
    }
    return kap;*/
    
    vector<real> gam;
    curvatures(m,kap,gam);
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

void ColocSim::evolve_system_RK4(real dp, bool fixdt) {
    // note, position_t is the heavy function here that solves the BEM problem.

    Mesh const& mesh = group.get_joined();

    Mesh m2 = mesh;
    Mesh m3 = mesh;
    Mesh m4 = mesh;
    
    CoordVec x1 = mesh.verts;
    PotVec p1 = get_phi();

    CoordVec k1_x = position_t(mesh,p1);

    real dt;
    if(fixdt) dt = dp;
    else dt = get_dt(dp,k1_x,kappa(mesh));

#ifdef VERBOSE
    cout << "\n dt = " << dt << endl << endl;
#endif
    CoordVec x2 = x1 + (0.5*dt)*k1_x;
    m2.verts = x2;
    PotVec k1_p = pot_t(mesh,k1_x);
    PotVec p2 = p1 + (0.5*dt)*k1_p;

    CoordVec k2_x = position_t(m2,p2);
    CoordVec x3 = x1 + (0.5*dt)*k2_x;
    m3.verts = x3;
    PotVec k2_p = pot_t(m2,k2_x);
    PotVec p3 = p1 + (0.5*dt)*k2_p;

    CoordVec k3_x = position_t(m3,p3);
    CoordVec x4 = x1 + dt*k3_x;
    m4.verts = x4;
    PotVec k3_p = pot_t(m3,k3_x);
    PotVec p4 = p1 + dt*k3_p;

    CoordVec k4_x = position_t(m4,p4);
    

    CoordVec xf = x1 + (dt/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
    

    //PotVec k4_p = pot_t(mesh,k4_x);
    //PotVec pf = p1 + (dt/6.0)*(k1_p + 2.0*k2_p + 2.0*k3_p + k4_p);

    PotVec pf = p1 + (dt/6.0)*(pot_t(mesh,k1_x) + 2.0*pot_t(m2,k2_x) + 2.0*pot_t(m3,k3_x) + pot_t(m4,k4_x));
    group.set_vertices(xf);
    set_phi(pf);
    time += dt;

    /*

    SimData x;
    x.pos = mesh.verts;
    x.phi = get_phi();
    
    SimData k1 = gradient(x);

    
    vector<real> kappa;
    generate_gradients(kappa);
    real dt = get_dt(dp,k1.pos,kappa);
    dt = 0.01;


    
    SimData k2 = gradient(x + 0.5*dt*k1);
    SimData k3 = gradient(x + 0.5*dt*k2);
    SimData k4 = gradient(x +     dt*k3);
    x += dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
    
    //x += dt*k1;


    mesh.verts = x.pos;
    set_phi(x.phi);
    time += dt;
    cout << "-----n\n\n    dt = " << dt << "\n\n-----u\n";
    */
}

void ColocSim::evolve_system(real dp, bool fixdt) {

    Eigen::MatrixXd G,H;
    assemble_matrices_prop(G,H);
    psi = solve_system(G,H*phi);

    //cout << setprecision(10) << psi << endl;

/*
    vector<vec3> vertex_gradients;
    vector<real> kappa;
    vertex_gradients = generate_gradients(kappa);

    vector<vec3> vertex_normals = mesh.generate_vertex_normals_max();

    for(size_t i(0);i<vertex_gradients.size();++i){
        vertex_gradients[i] = vertex_normals[i]*psi(i);
    }*/

    Mesh mesh = group.get_joined();

    vector<vec3> normals = generate_triangle_normals(mesh); // to ensure there are normals
    vector<vector<size_t>> triangle_indices = generate_triangle_indices(mesh);
    vector<vec3> triangle_gradients = generate_tangent_gradients(mesh,make_copy(phi));

    vector<vec3> vertex_normals = generate_vertex_normals(mesh);

    vector<vec3> vertex_gradients;


#if LINEAR
    for(size_t i(0);i<mesh.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {

            grad += triangle_gradients[index]+psi(i)*normals[index];
            //grad += triangle_gradients[index];
            num++;
        }
        grad *= 1.0/num;

        //grad -= grad.dot(vertex_normals[i])*vertex_normals[i]; // remove normal part
        //grad += vertex_normals[i]*psi(i);

        vertex_gradients.push_back(grad);

        //cout << grad.norm() << endl;
    }
#else
    vertex_gradients = vector<vec3>(mesh.verts.size());
    vector<real> num(mesh.verts.size(),0.0);
    for(size_t i(0);i<mesh.trigs.size();++i) {
        Triplet t(mesh.trigs[i]);
        Cubic interp(mesh.verts[t.a],mesh.verts[t.b],mesh.verts[t.c],vertex_normals[t.a],vertex_normals[t.b],vertex_normals[t.c]);
        num[t.a]++; num[t.b]++; num[t.c]++;
        vertex_gradients[t.a] += interp.tangent_derivative_at_a(phi(t.a),phi(t.b),phi(t.c));
        vertex_gradients[t.b] += interp.tangent_derivative_at_b(phi(t.a),phi(t.b),phi(t.c));
        vertex_gradients[t.c] += interp.tangent_derivative_at_c(phi(t.a),phi(t.b),phi(t.c));
    }
    for(size_t i(0);i<mesh.verts.size();++i) {
        vertex_gradients[i] = vertex_gradients[i]*(1.0/num[i]) + psi(i)*vertex_normals[i];
    }
#endif

    vector<real> kappa, gamma;
    //generate_curvature_tensor(mesh,kappa,gamma);
    //vertex_gradients = 
    generate_gradients(kappa);

    real dt;
    if(fixdt) dt = dp;
    else dt = get_dt(dp,vertex_gradients,kappa);
    
#ifdef VERBOSE
    cout << "\n dt = " << dt << endl << endl;
#endif
    
    //cout << psi << endl;
    
    for(size_t i(0);i<mesh.verts.size();++i) {
        mesh.verts[i] += vertex_gradients[i]*dt;

        //vec3 norm(x[i]);
        //norm.normalize();
        //cout << vertex_normals[i].dot(norm) << endl;
    }
    
    group.set_vertices(mesh.verts);

    //mesh.relax_vertices();

    // recompute kappa and vertex grads after volume change
    // vertex_gradients = generate_gradients(kappa);
    //generate_curvature_tensor(mesh,kappa,gamma);
    //vertex_gradients = 
    generate_gradients(kappa);

    real vol = volume(mesh);
    for(size_t i(0);i<mesh.verts.size();++i) {
        phi(i) += dt*potential_t(vertex_gradients[i].norm2(),vol,kappa[i]);

        //psi(i) = kappa[i];
    }
    
    
    time += dt;

}



vector<vec3> ColocSim::generate_gradients(vector<real>& kappa,Mesh const& m,Eigen::VectorXd const& phi,Eigen::VectorXd const& psi) const {

    kappa.clear();
    vector<vec3> gradients;
    const vector<vec3>& x = m.verts;
    vector<vector<size_t>> verts = generate_neighbours(m);

    for(size_t i(0);i<verts.size();++i) {
        verts[i].push_back(i);
    }
    
    vector<vec3> norms = generate_vertex_normals(m);

    size_t n(x.size());

    for(size_t i(0);i<n;++i) {
        vector<vec3> positions;
        vector<vec3> potentials;
        for(size_t j : verts[i]) {
            vec3 pos(x[j]);
            positions.push_back(pos);
            pos.z = phi(j);
            potentials.push_back(pos);
        }
        FittingTool fit;
        fit.compute_quadratic_fit(norms[i],x[i],positions);
        vec3 normal(fit.get_normal());
        kappa.push_back(fit.get_curvature());

        // I think this is not so rigorous (see Wang 2014)
        fit.compute_quadratic_fit(normal,x[i],potentials);
        vector<real> params(fit.get_params());
        CoordSystem system(fit.copy_coord_system());
        vec3 tangrad = system.world_coords_relative(params[1],params[2],0.0);

        gradients.push_back(tangrad + psi(i)*normal);
        
    }

    return gradients;
}


vector<vec3> ColocSim::generate_gradients(vector<real>& kappa) {

    
    
#if ANALYTICS
    Mesh output_mesh_coordsys;
    Mesh output_mesh_approx;
    real factor(0.05);
#endif
    
    Mesh const& mesh = group.get_joined();

    kappa.clear();
    vector<vec3> gradients;
    const vector<vec3>& x = mesh.verts;
    vector<vector<size_t>> verts = generate_2_ring(mesh);
    // if only 1_ring:
    
    
    verts = generate_neighbours(mesh);

    for(size_t i(0);i<verts.size();++i) {
        verts[i].push_back(i);
    }
    


    vector<vector<size_t>> nbs = generate_neighbours(mesh);
    vector<real> d(mesh.verts.size(),0.0);
    for(size_t i(0);i<mesh.verts.size();++i) {
        for(size_t j : nbs[i]) {
            d[i] += (x[j]-x[i]).norm();
        }
        d[i] /= static_cast<real>(nbs[i].size());
        //cout << d[i] << endl;
    }
    
    
    vector<vec3> norms = generate_vertex_normals(mesh);

    size_t n(x.size());

    for(size_t i(0);i<n;++i) {
        vector<vec3> positions;
        vector<vec3> potentials;
        for(size_t j : verts[i]) {
            vec3 pos(x[j]);
            positions.push_back(pos);
            pos.z = phi(j);
            potentials.push_back(pos);
        }
        FittingTool fit;
        fit.compute_quadratic_fit(norms[i],x[i],positions);
        vec3 normal(fit.get_normal());
        kappa.push_back(fit.get_curvature());

        //cout << fit.get_curvature() << ' ';

        // I think this is not so rigorous (see Wang 2014)
        fit.compute_quadratic_fit(normal,x[i],potentials);
        vector<real> params(fit.get_params());
        CoordSystem system(fit.copy_coord_system());
        vec3 tangrad = system.world_coords_relative(params[1],params[2],0.0);

        gradients.push_back(tangrad + psi(i)*normal);
        //cout << " tan: " << tangrad << " + norm: " << psi(i)*normal << "  psi = " << psi(i) << " normal = " << normal << endl;

#if ANALYTICS
        Mesh coordsys;
        coordsys.add_vertex(system.world_coords(0.0,0.0,0.0));
        coordsys.add_vertex(system.world_coords(factor,0.0,0.0));
        coordsys.add_vertex(system.world_coords(0.0,factor,0.0));
        coordsys.add_vertex(system.world_coords(0.0,0.0,factor));

        coordsys.add_triangle(0,1,2);
        coordsys.add_triangle(0,2,3);
        coordsys.add_triangle(0,3,1);

        Mesh approx;

        float dx = 0.03;
        float extent = 0.3;
        size_t ns = 0,nt = 0;
        for(float s(-extent);s<=extent;s+=dx) {
            ns++;
            for(float t(-extent);t<=extent;t+=dx) {
                approx.add_vertex(fit.get_position(s,t));
            }
        }
        nt = approx.verts.size()/ns;

        for(size_t i(1);i<ns;++i) {
            for(size_t j(1);j<nt;++j) {
                approx.add_triangle((i-1)*nt + (j-1),i*nt + (j-1),i*nt + j);
                approx.add_triangle((i-1)*nt + (j-1),i*nt + j,(i-1)*nt + j);
            }
        }

        output_mesh_coordsys.add_mesh(coordsys);
        output_mesh_approx.add_mesh(approx);
#endif
        
    }

#if ANALYTICS
    output_mesh_coordsys.export_ply("test-coord-sys.ply");
    output_mesh_approx.export_ply("test-approx.ply");
#endif
    return gradients;
}


vector<vec3> ColocSim::generate_tangent_gradients(Mesh const& m, vector<real> const& pot) const {
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

vector<real> ColocSim::curvature_param() {
    Mesh const& mesh = group.get_joined();

    vector<real> max_curv = max_curvature(mesh);

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
        max_curv[i] += 1.0*meandiffgrad[i];
    }

    // smoothing the parameters
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

void ColocSim::remesh(real L) {
    
}

void ColocSim::remesh(Mesh& mesh,vector<real>& pot,real L) {
    
    
    PotVec new_curv_params = curvature_param();


    curvature_params = damping_factor*curvature_params + (1.0-damping_factor)*new_curv_params;

    
    HalfedgeMesh manip(generate_halfedges(mesh));

    split_edges(manip,curvature_params,L*3.0/4.0);
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
    collapse_edges(manip,curvature_params,L*4.0/5.0);
    flip_edges(manip,1);
    flip_edges(manip,1);
    relax_vertices(manip);
    flip_edges(manip,1);
    relax_vertices(manip);
    vector<real> new_pot;
    Mesh new_mesh = generate_mesh(manip);
    project_and_interpolate(new_mesh,new_pot, mesh, pot);
    mesh = new_mesh;
    pot = new_pot
    
}



void ColocSim::export_mesh_psi(string fname) const {
    export_ply_float(fname,mesh,make_copy(psi));
}

void ColocSim::export_mesh(string fname) const {
    vector<real> col;
    Eigen::VectorXd vec = phi;
    real min_value(vec(0));
    real max_value(vec(0));
    for(size_t i(0);i<mesh.verts.size();++i) {
        col.push_back(vec(i));
        min_value = min(min_value,vec(i));
        max_value = max(max_value,vec(i));
    }
#ifdef VERBOSE
    cout << "interval: [" << min_value << ',' << max_value << "]." << endl;
#endif
    export_ply_float(fname,mesh,col);
}


} // namespace Group

} // namespace Bem