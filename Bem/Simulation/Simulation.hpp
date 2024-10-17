#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <string>
#include <vector>

#include "../basic/Bem.hpp"
#include "../Mesh/Mesh.hpp"
#include "../Integration/Integrator.hpp"

#include <Eigen/Dense>

#include <cassert>

namespace Bem {


// functions for translating between PotVec and Eigen::VectorXd
std::vector<real> make_copy(Eigen::VectorXd const& vec);
Eigen::VectorXd make_copy(std::vector<real> const& vec);

// default pressure field
real default_field(vec3 x,real t);

class Simulation {
public:

    Simulation(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :p_inf(p_inf), 
        epsilon(epsilon),  
        sigma(sigma), 
        gamma(gamma),
        pressurefield(pressurefield),
        time(0.0),
        min_dt(-1.0),
        dp_balance(3.0),
        bicgstab(true),
        num_threads(100),
        mesh(initial) {
            // initializing other default values:
            V_0 = volume(initial);

            inter.set_quadrature(quadrature_7);
            inter.set_quadrature(gauss_7);

#ifdef VERBOSE
            std::cout << "Simulation initialized with N = " << mesh.verts.size() << " vertices and M = " << mesh.trigs.size() << " triangles." << std::endl;
            std::cout << " p_inf:   " << p_inf << std::endl;
            std::cout << " epsilon: " << epsilon << std::endl;
            std::cout << " sigma:   " << sigma << std::endl;
            std::cout << " gamma:   " << gamma << std::endl;
#endif
        }

    virtual ~Simulation() {}

    // dimensions of the approximation spaces describing phi resp. psi
    virtual size_t phi_dim() const = 0;
    virtual size_t psi_dim() const = 0;

    // function that fills the system matrices G and H corresponding to a mesh m.
    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const = 0;
    // if no other mesh is specified, the mesh describing the current state of the simulation is used:
    void assemble_matrices_prop(Eigen::MatrixXd& G,Eigen::MatrixXd& H) const {
        assemble_matrices(G,H,mesh);
    }

    // solving the system G*psi = H_phi = H*phi for psi (returned vector)
    Eigen::VectorXd solve_system(Eigen::MatrixXd const& G,Eigen::VectorXd const& H_phi) const;

    // This function only computes the psi values without evolving the system in time
    void compute_psi() {
        Eigen::MatrixXd G,H;
        assemble_matrices_prop(G,H);
        psi = solve_system(G,H*phi);
    }

    // function that handels the time evolution of the system: must be defined in subclasses
    virtual void evolve_system(real dp, bool fixdt = false) = 0;


    // function to export the mesh representing the bubble's current geometry together
    // with the values of phi and psi on its vertices
    void export_mesh(std::string fname) const;
    // this function allows to export another set of scalar vertex data instead of phi and psi
    void export_mesh_values(std::string fname,std::vector<real> values) const;


    // getter and setter functions to get/set the main state variables of the system:

    void set_phi(std::vector<real> const& value) {
        assert(value.size() == phi_dim());
        phi = make_copy(value);
    }

    void set_psi(std::vector<real> const& value) {
        assert(value.size() == psi_dim());
        psi = make_copy(value);
    }

    void set_phi(real value) {
        phi.fill(value);
    } 

    void set_V_0(real value) {
        V_0 = value;
    }

    real get_phi(size_t index) const {
        assert(index < phi_dim());
        return phi(index);
    }

    real get_psi(size_t index) const {
        assert(index < psi_dim());
        return psi(index);
    }

    void set_time(real const& value) {
        time = value;
    }
    
    real get_time() const {
        return time;
    }

    real get_volume() const {
        return volume(mesh);
    }

    
    const std::vector<vec3>& get_vertices() {
        return mesh.verts;
    }

    std::vector<real> get_psi() {
        return make_copy(psi);
    }

    std::vector<real> get_phi() {
        return make_copy(phi);
    }

    void set_min_dt(real value) {
        min_dt = value;
    }

    void set_dp_balance(real value) {
        dp_balance = value;
    }

    void set_bcgstab(bool value) {
        bicgstab = value;
    }

    void set_num_threads(size_t num) {
        num_threads = num;
    }

    void set_quadrature(std::vector<quadrature_2d> const& quad) {
        inter.set_quadrature(quad);
    }

    void set_quadrature(std::vector<quadrature_1d> const& quad) {
        inter.set_quadrature(quad);
    }

protected:

    real get_dt(real dp,std::vector<vec3> const& gradients, std::vector<real> const& grad_potential) const;
    // time derivative of the potential
    real potential_t(real grad_squared, real volume, real kappa, vec3 pos, real t) const;


    // The following constants are given in simulation units;
    // Eventual conversions happen outside this class

    // necessary constants (to run the simulation)
    real p_inf, epsilon, sigma, gamma;

    // initial condition informations:
    real V_0;

    // wave informations:
    real (*pressurefield)(vec3 x,real t);

    // simulation time:
    real time;
    real min_dt;
    real dp_balance;

    // which solver
    bool bicgstab;

    // number of threads for matrix generation (if supported)
    size_t num_threads;

    // Mesh describing the Bubble surface and an integrator object,
    // providing some functions to integrate over the mesh.
public:
    Mesh mesh; 
protected:
    Integrator inter;

    // Vectors that approximate the potential (phi) and its normal derivative (psi) on the surface
    Eigen::VectorXd psi;
    Eigen::VectorXd phi;

};

} // namespace Bem

#endif // SIMULATION_HPP