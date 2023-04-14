#include "Simulation.hpp"
#include "../Integration/Integrator.hpp"
#include "../Integration/ResultTypes.hpp"
#include "../basic/Bem.hpp"
#include "../Mesh/HalfedgeMesh.hpp"
#include "../Mesh/MeshManip.hpp"
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

real default_waveform(vec3 x,real t) {
    return x.dot(vec3()) + t*0.0; // avoid warnings for not using x and t
}

Eigen::VectorXd Simulation::solve_system(Eigen::MatrixXd const& G,Eigen::VectorXd const& H_phi) const {
    // attention! right now, G is completely symmetric! could only use one half!
#ifdef VERBOSE
    cout << " solving system..." << flush;
    auto start = high_resolution_clock::now();
#endif
    // PartialPivLU needs an invertible SQUARE matrix! to be sure, can use FullPivLU, but not in parallel!

    Eigen::VectorXd x;
    if(bicgstab){
        Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
        solver.compute(G);
        x = solver.solve(H_phi);
        //Eigen::VectorXd x = solver.solveWithGuess(H_phi,psi); // leads to instabilities in computation...
    } else {
        Eigen::PartialPivLU<Eigen::MatrixXd> solver;
        solver.compute(G);
        x = solver.solve(H_phi);
    }
    
/*
    Eigen::ConjugateGradient<Eigen::MatrixXd,Eigen::Lower|Eigen::Upper> cg; 
    cg.compute(G);
    Eigen::VectorXd x = cg.solve(H_phi);
    */

#ifdef VERBOSE
    cout << " - done." << endl;
    auto end = high_resolution_clock::now();
    cout << "used time = " << duration_cast<duration<double>>(end-start).count() << " s. " << endl;
    //cout << "estimated error: " << solver.error() << endl;
    cout << "true error: " << (G*x - H_phi).maxCoeff() << endl;
#endif
    return x;
}

real Simulation::potential_t(real grad_squared, real volume, real kappa,vec3 pos, real t) const {
    /*
    cout << "---- vals ----" << endl;
    cout << 2.0*sigma*kappa << endl;
    cout << 0.5*grad_squared << endl;
    cout << p_inf << endl;
    cout << epsilon*pow(V_0/volume,gamma) << endl;
    cout << waveform(pos,t) << endl;
    cout << 2.0*sigma*kappa + 0.5*grad_squared + p_inf - epsilon*pow(V_0/volume,gamma) + waveform(pos,t) << endl;
    */
    return  2.0*sigma*kappa + 0.5*grad_squared + p_inf - epsilon*pow(V_0/volume,gamma) + waveform(pos,t);
}

real Simulation::get_dt(real dp,std::vector<vec3> const& gradients, std::vector<real> const& grad_potential) const {
    real max_pot(0.0);
    real max_vel(0.0);
    for(size_t i(0);i<gradients.size();++i) {
        max_pot = max(abs(grad_potential[i]),max_pot);
        max_vel = max(gradients[i].norm2(),max_vel);
    }
    max_vel = sqrt(max_vel);

    // dp_balance is parameter to balance the influence of the time- and space derivatives on the timestep
    real dt = dp/(max_pot + dp_balance*max_vel);

    if(min_dt > 0.0) dt = min(dt,min_dt); // if min_dt is specified, take the minimum!

    return dt;
}

void Simulation::export_mesh(string fname) const {
    export_ply_float(fname,mesh,make_copy(phi),make_copy(psi));
}

void Simulation::export_mesh_values(string fname,vector<real> values) const {
    export_ply_float(fname,mesh,values);
}



} // namespace Bem