#include "ConLinSim.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
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



void ConLinSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
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


    return;

    /*
    
    const vector<vec3>& x(m.verts);


    for(size_t i(0);i<m.trigs.size();++i) {
        for(size_t j(0);j<m.trigs.size();++j) {
            inter.integrate_ConLin(x,m.trigs[i],m.trigs[j],i,j,G,H);
        }
#ifdef VERBOSE
        cout << " Assembling matrices... progress: " << float(i+1)/m.trigs.size()*100.0 << "%                                    \r" << flush;
#endif
    }
#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif

    */
}

CoordVec ConLinSim::position_t(Mesh const& m,PotVec const& pot) const {

    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);
    vector<vec3> normals = generate_triangle_normals(m);
    vector<vec3> tangent_gradients = generate_tangent_gradients(m,pot);
    vector<vec3> vertex_gradients;

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

vector<vec3> ConLinSim::generate_tangent_gradients(Mesh const& m, vector<real> const& pot) const {
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

/*
void ConLinSim::evolve_system(real dp) {

    Eigen::MatrixXd G,H;
    assemble_matrices_prop(G,H);
    psi = solve_system(G,H*phi);
    
    vector<vec3> normals = mesh.generate_triangle_normals(); // to ensure there are normals
    vector<vector<size_t>> triangle_indices = mesh.generate_triangle_indices();
    vector<vec3> triangle_gradients = generate_tangent_gradients();
    const vector<vec3>& x(mesh.get_vertices());

    vector<vec3> vertex_gradients;
    for(size_t i(0);i<mesh.N();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {
            // '-' because psi is solution of problem with inward pointing normals
            grad += triangle_gradients[index]-psi(index)*normals[index];
            num++;
        }
        grad *= 1.0/num;

        vertex_gradients.push_back(grad);
    }

    real dt = get_dt(dp,vertex_gradients);

    // update vertex positions:
    Eigen::VectorXd vec = phi;
    real min_value(vec(0));
    real max_value(vec(0));
    for(size_t i(0);i<mesh.N();++i) {
        mesh.set_vertex_position(x[i]+vertex_gradients[i]*dt,i);
        min_value = min(min_value,vec(i));
        max_value = max(max_value,vec(i));
        phi(i) += dt*potential_t(vertex_gradients[i].norm2());
    }
    time += dt;
#ifdef VERBOSE
    cout << " - range = [" << min_value << ',' << max_value << ']' << endl;
#endif
    // update potential values:
    //evolve_potential(dt,vertex_gradients);
}



vector<vec3> ConLinSim::generate_tangent_gradients() const {
    vector<vec3> gradients;
    const vector<vec3>& x = mesh.get_vertices();
    const vector<Triplet>& trigs = mesh.get_triangles();
    for(Triplet trig : trigs) {
        // this formula for the gradient of linear functions on 
        // an unstructured grid is taken from: Dombre_2019, page 8.
        vec3 ab(x[trig.b]-x[trig.a]);
        vec3 bc(x[trig.c]-x[trig.b]);
        
        vec3 n(ab.vec(bc));
        n *= 1.0/n.norm2();

        vec3 grad = (phi(trig.c)-phi(trig.b))*n.vec(ab)
                  + (phi(trig.a)-phi(trig.b))*n.vec(bc);

        gradients.push_back(grad);
    }
    return gradients;
}

void ConLinSim::export_mesh_psi(string fname) const {
    std::vector<real> values;
    for(size_t i(0);i<mesh.M();++i){
        for(size_t j(0);j<3;++j){
            values.push_back(psi(i));
        }
    }
    mesh.export_ply_float_separat(fname,values);
}

void ConLinSim::export_mesh(string fname) const {
    vector<real> col;
    Eigen::VectorXd vec = phi;
    real min_value(vec(0));
    real max_value(vec(0));
    for(size_t i(0);i<mesh.N();++i) {
        col.push_back(vec(i));
        min_value = min(min_value,vec(i));
        max_value = max(max_value,vec(i));
    }
#ifdef VERBOSE
    cout << "interval: [" << min_value << ',' << max_value << "]." << endl;
#endif
    mesh.export_ply_float(fname,col);
}
*/