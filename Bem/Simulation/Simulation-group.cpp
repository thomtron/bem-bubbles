#include "Simulation-group.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
#include "../basic/Bem.hpp"
#include "../Mesh/MeshIO.hpp"
#include <vector>
#ifdef VERBOSE
#include <chrono>
#endif

#include <Eigen/IterativeLinearSolvers> // for cg solver

using namespace std;
#ifdef VERBOSE
using namespace chrono;
#endif


namespace Bem {

namespace Group {

vector<real> make_copy(Eigen::VectorXd const& vec) {
    return vector<real>(vec.begin(),vec.end());
}

Eigen::VectorXd make_copy(vector<real> const& vec) {
    size_t N(vec.size());
    Eigen::VectorXd result(N);
    for(size_t i(0);i<N;++i){
        result(i) = vec[i];
    }
    return result;
}

Eigen::VectorXd Simulation::solve_system(Eigen::MatrixXd const& G,Eigen::VectorXd const& H_phi) const {
    // attention! right now, G is completely symmetric! could only use one half!
#ifdef VERBOSE
    cout << " solving system..." << flush;
    auto start = high_resolution_clock::now();
#endif
    // PartialPivLU needs an invertible SQUARE matrix! to be sure, can use FullPivLU, but not in parallel!
    Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
    solver.compute(G);
    //Eigen::VectorXd x = solver.solveWithGuess(H_phi,psi);
    //Eigen::FullPivLU<Eigen::MatrixXd> solver;
    //solver.compute(G);
    Eigen::VectorXd x = solver.solve(H_phi);
    
/*
    Eigen::ConjugateGradient<Eigen::MatrixXd,Eigen::Lower|Eigen::Upper> cg; 
    cg.compute(G);
    Eigen::VectorXd x = cg.solve(H_phi);
    */

    //cout << H*phi -G*psi<< endl;
#ifdef VERBOSE
    cout << " - done." << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
    //cout << "estimated error: " << solver.error() << endl;
    cout << "true error: " << (G*x - H_phi).maxCoeff() << endl;
#endif
    return x;
}

real Simulation::potential_t(real grad_squared) const {
    // kappa term not implemented here -> depends on local 
    // geometry -> additional argument necessary! (index of element etc.)

    // taib ?!
    //real kappa = 0.0;
    //return p_inf - (epsilon - 2.0*sigma*kappa - (1.0 - epsilon + 2.0*sigma*kappa_0)*pow(V_0/mesh.volume(),gamma)) + 0.5*grad_squared;


    // Wang_2014
    //return 1.0 - epsilon*pow(V_0/mesh.volume(),gamma)  + 0.5*grad_squared;

    // we take here p_vap = epsilon, rho = 1.0
    return p_inf - epsilon + 0.5*grad_squared;
    //return 0.3 + 0.5*grad_squared-1.0/mesh.volume();
}

real Simulation::potential_t(real grad_squared,real kappa) const {
    real V0(0.0);
    real V(0.0);
    for(real elm : V_0)
        V0 += elm;
    for(real elm : group.volumes())
        V += elm;
    return potential_t(grad_squared,V/V0,kappa);
}

real Simulation::potential_t(real grad_squared, real volume, real kappa) const {
    // Wang_2014 II
    return  2.0*sigma*kappa + 0.5*grad_squared + 1.0 - epsilon*pow(1.0/volume,lambda);
}

real Simulation::potential_t(real grad_squared, real volume, real kappa, real t) const {
    // Wang_2014 II
    return  2.0*sigma*kappa + 0.5*grad_squared + 1.0 - epsilon*pow(1.0/volume,lambda) + 0.5*sin(t);
}


real Simulation::get_dt(real dp,vector<vec3> const& gradients) const {
    real max_value(0.0);
    real max_velocity(0.0);
    for(size_t i(0);i<phi_dim();++i) {
        max_value = max(max_value,abs(potential_t(gradients[i].norm2())));
        max_velocity = max(max_velocity,gradients[i].norm());
    }
    return dp/(max_value+max_velocity);
}

real Simulation::get_dt(real dp,vector<vec3> const& gradients, vector<real> const& kappa) const {
    real max_value(0.0);
    real max_velocity(0.0);
    vector<real> vols = group.volumes();
    for(size_t i(0);i<vols.size();++i)
        vols[i] /= V_0[i];
    vols = group.expand_to_vertex_data(vols);

    for(size_t i(0);i<phi_dim();++i) {
        max_value = max(max_value,abs(potential_t(gradients[i].norm2(),vols[i],kappa[i])));
        max_velocity = max(max_velocity,gradients[i].norm());
    }
    return dp/(max_value+max_velocity);
}

real Simulation::get_dt(real dp,vector<vec3> const& gradients, vector<real> const& kappa, real t) const {
    real max_value(0.0);
    real max_velocity(0.0);
    vector<real> vols = group.volumes();
    for(size_t i(0);i<vols.size();++i)
        vols[i] /= V_0[i];
    vols = group.expand_to_vertex_data(vols);
    
    for(size_t i(0);i<phi_dim();++i) {
        max_value = max(max_value,abs(potential_t(gradients[i].norm2(),vols[i],kappa[i],t)));
        max_velocity = max(max_velocity,gradients[i].norm());
    }
    return dp/(max_value+max_velocity);
}

void Simulation::export_mesh(string fname) {
    export_ply(fname,group.get_joined());
}

void Simulation::export_mesh_values(string fname,vector<real> values) {
    export_ply_float(fname,group.get_joined(),values);
}

} // namespace Group

} // namespace Bem