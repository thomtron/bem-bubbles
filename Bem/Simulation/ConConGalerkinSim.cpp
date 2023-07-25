#include "ConConGalerkinSim.hpp"
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


void ConConGalerkinSim::assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const {
#ifdef VERBOSE
    auto start = high_resolution_clock::now();
#endif
    G = Eigen::MatrixXd::Zero(m.trigs.size(),m.trigs.size());
    H = Eigen::MatrixXd::Zero(m.trigs.size(),m.trigs.size());

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
    
    Mesh local(m);
    const CoordVec& x(local.verts);
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

#ifdef VERBOSE
    cout << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
#endif

    return;
}


CoordVec ConConGalerkinSim::position_t(Mesh const& m,PotVec const& pot) const {
    
    Eigen::MatrixXd G,H;
    assemble_matrices(G,H,m);
    Eigen::VectorXd psi_l = solve_system(G,H*make_copy(pot));

    vector<vector<size_t>> triangle_indices = generate_triangle_indices(m);
    CoordVec normals = generate_triangle_normals(m);
    CoordVec tangent_gradients = generate_tangent_gradients(m,pot);
    CoordVec vertex_gradients;

    // averaging the normal derivatives to obtain their value at the vertices and then 
    // adding the tangent derivatives of phi to the normal derivatives.
    for(size_t i(0);i<m.verts.size();++i) {
        real num(0.0);
        vec3 grad;
        for(size_t index : triangle_indices[i]) {

            grad += psi_l(index)*normals[index];
            num++;
        }
        grad *= 1.0/num;
        grad += tangent_gradients[i];
        // note a potential problem with this code: the tangent gradients are not necessarily
        // orthogonal to grad, as grad is composed of the normals of the adjacent triangles and 
        // does not in general point into the same direction as the normal direction obtained by
        // the surface fit. An alternative would be to store the vertex normal from the surface 
        // fit, averaging the psi_l values at the vertex and then multiplying this average with 
        // the vertex normal.

        vertex_gradients.push_back(grad);
    }

    return vertex_gradients;
}

#include <set>

// compute approximation to per vertex tangent gradient using a local fit of the mesh
CoordVec ConConGalerkinSim::generate_tangent_gradients(Mesh const& m, PotVec const& pot) const {
    
    vector<vector<size_t>> neighbours = generate_neighbours(m);
    vector<vector<size_t>> trig_inds = generate_triangle_indices(m);

    vector<vector<size_t>> verts; // to include second neighbouring triangles 
                                  // otherwise, just take verts = generate_triangle_indices(m)
                                  // Note, that we need triangle indices here and thus cannot
                                  // simply use the function generate_2_ring()
    for(vector<size_t> const& elm : neighbours) {
        set<size_t> inds;
        for(size_t k : elm) {
            for(size_t j : trig_inds[k]) {
                inds.insert(j);
            }
        }
        verts.push_back(vector<size_t>(inds.begin(),inds.end()));
    }


    CoordVec norms = generate_vertex_normals(m);


    CoordVec gradients;
    for(size_t i(0);i<m.verts.size();++i) {

        // Here you have to ensure that enopugh points are used! I've
        // added the point at the vertex with the mean potential of the 
        // triangles around it, but in general there may be needed more!
        CoordVec positions;
        PotVec potentials;
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
        // to ensure that the normal corresponds exactly to the vertical axis:
        fit.compute_quadratic_fit(normal,m.verts[i],positions); 
        
        // transforming the positions into the local coordinate system and replacing
        // their z-values with the phi-values for obtaining below an approximate fit 
        // of phi on the surface. Note that this approximation is only meaningful when 
        // the surface is relatively flat.
        CoordSystem system(fit.copy_coord_system());
        for(size_t j(0);j<positions.size();++j) {
            vec3 trans = system.transform(positions[j]);
            trans.z = potentials[j];
            positions[j] = trans;
        }

        // fit of the potential values
        fit.compute_quadratic_fit(vec3(0.0,0.0,1.0),vec3(),positions);
        vector<real> params(fit.get_params());
        CoordSystem system2(fit.copy_coord_system());

        // pushing the gradient of phi(x,y) at (0,0,0): (params[1],params[2],0.0) back to the 
        // world coordinate system
        vec3 tangrad = system.world_coords_relative(system2.world_coords_relative(params[1],params[2],0.0));

        gradients.push_back(tangrad);
    }
    return gradients;
}
