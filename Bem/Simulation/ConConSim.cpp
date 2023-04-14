#include "ConConSim.hpp"
#include "../Integration/Integrator.hpp"
#include "../Mesh/FittingTool.hpp"
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

// NOTE!! HIER SIND IMMER NOCH VERALTETE FUNKTIONEN ENTHALTEN - UPDATEN BEVOR INS CMAKELISTS INCLUDEN!


void ConConSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif
    G = Eigen::MatrixXd::Zero(m.trigs.size(),m.trigs.size());
    H = Eigen::MatrixXd::Zero(m.trigs.size(),m.trigs.size());

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
        
        for(size_t i(0);i<M;++i) {
            int_local.integrate_Con(x,local.trigs[i],local.trigs[j],G(i,j),H(i,j));
        }

#ifdef VERBOSE
        if(omp_get_thread_num() == 0)
            cout << " Assembling matrices... progress (approx.): " << float(j+1)/M*100.0*omp_get_num_threads() << "%                                    \r" << flush;
#endif 
        
    }

    }


    return;





    /*

    for(size_t i(0);i<m.trigs.size();++i) {
        for(size_t j(0);j<m.trigs.size();++j) {
            inter.integrate_Con(m.verts,m.trigs[i],m.trigs[j],G(i,j),H(i,j));
        }
#ifdef VERBOSE
        cout << " Assembling matrices... progress: " << float(i+1)/m.trigs.size()*100.0 << "%                                    \r" << flush;
#endif
    }
#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif*/
}


CoordVec ConConSim::position_t(Mesh const& m,PotVec const& pot) const {
    
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

            grad += psi_l(index)*normals[index];
            num++;
        }
        grad *= 1.0/num;
        grad += tangent_gradients[i];

        vertex_gradients.push_back(grad);
    }

    return vertex_gradients;

}

#include <set>

// compute approximation to per vertex tangent gradient
vector<vec3> ConConSim::generate_tangent_gradients(Mesh const& m, vector<real> const& pot) const {
    
    vector<vector<size_t>> neighbours = generate_neighbours(m);
    vector<vector<size_t>> trig_inds = generate_triangle_indices(m);

    vector<vector<size_t>> verts; // to include second neighbouring triangles 
                                  // otherwise, just take verts = generate_triangle_indices(m)
    for(vector<size_t> const& elm : neighbours) {
        set<size_t> inds;
        for(size_t k : elm) {
            for(size_t j : trig_inds[k]) {
                inds.insert(j);
            }
        }
        verts.push_back(vector<size_t>(inds.begin(),inds.end()));
    }


    vector<vec3> norms = generate_vertex_normals(m);


    vector<vec3> gradients;
    for(size_t i(0);i<m.verts.size();++i) {

        // Here you have to ensure that enopugh points are used! I've
        // added the point at the vertex with the mean potential of the 
        // triangles around it, but in general there may be needed more!
        vector<vec3> positions;
        vector<real> potentials;
        real mean_pot(0.0);
        for(size_t j : verts[i]) {
            Triplet t(m.trigs[j]);
            vec3 pos = (1.0/3.0)*(m.verts[t.a]+m.verts[t.b]+m.verts[t.c]);
            positions.push_back(pos);
            mean_pot += pot[j];
            potentials.push_back(pot[j]);
        }
        mean_pot /= static_cast<real>(verts[i].size());
        positions.push_back(m.verts[i]);
        potentials.push_back(mean_pot);

        FittingTool fit;
        fit.compute_quadratic_fit(norms[i],m.verts[i],positions);
        vec3 normal(fit.get_normal());
        fit.compute_quadratic_fit(normal,m.verts[i],positions); // just to have exactly the same coord system as below
        
        //cout << endl << '(';
        CoordSystem system(fit.copy_coord_system());
        for(size_t j(0);j<positions.size();++j) {
            vec3 trans = system.transform(positions[j]);
            trans.z = potentials[j];
            positions[j] = trans;

            //cout << trans << ") (";
        }
        //cout << ')' << endl;


        // I think this is not so rigorous (see Wang 2014)
        fit.compute_quadratic_fit(vec3(0.0,0.0,1.0),vec3(),positions); // now with potential values
        vector<real> params(fit.get_params());
        CoordSystem system2(fit.copy_coord_system());

        //for(real elm : params)
        //    cout << elm << " _ ";
        //cout << endl;
        vec3 tangrad = system.world_coords_relative(system2.world_coords_relative(params[1],params[2],0.0));

        //cout << tangrad << endl;

        gradients.push_back(tangrad);
    }
    return gradients;
}


/*
void ConConSim::evolve_system(real dp) {

    real dt = dp; // Note, here one could insert adaptive timestepping!

    vector<vec3> gradients;

    mesh.generate_triangle_normals(); // to ensure there are normals
    for(size_t i(0);i<mesh.N();++i) {
        vec3 approx_gradient;
        float num(0.0);
        for(size_t t : mesh.triangle_indices[i]){
            // "-=" due to normals that point to exteriour
            approx_gradient -= mesh.triangle_normals[t]*psi(t);
            num++;
        }
        approx_gradient *= (1.0/num);

        gradients.push_back(approx_gradient);

        mesh.vertices[i] += approx_gradient*dt;
    }
    //evolve_potential(dt,gradients);
}
*/
