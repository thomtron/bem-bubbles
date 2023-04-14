#ifndef SIMULATION_GROUP_HPP
#define SIMULATION_GROUP_HPP

#include <iostream>
#include <string>
#include <vector>

#include "../basic/Bem.hpp"
#include "../Mesh/Mesh.hpp"
#include "../Integration/Integrator.hpp"
#include "SimData.hpp"
#include "MeshGroup.hpp"

#include <Eigen/Dense>

#include <cassert>

namespace Bem {

namespace Group {

template<typename func_t,typename array_t>
array_t RK4(func_t f, array_t const& x,real t,real dt) {
    array_t k1 = f(x,t);
    array_t k2 = f(x + 0.5*dt*k1,t + 0.5*dt);
    array_t k3 = f(x + 0.5*dt*k2,t + 0.5*dt);
    array_t k4 = f(x +     dt*k3,t +     dt);
    return x + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
}

template<typename func_t,typename array_t>
array_t RK4(func_t f, array_t const& x, real dt) {
    array_t k1 = f(x);
    array_t k2 = f(x + 0.5*dt*k1);
    array_t k3 = f(x + 0.5*dt*k2);
    array_t k4 = f(x +     dt*k3);
    return x + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
}

std::vector<real> make_copy(Eigen::VectorXd const& vec);

Eigen::VectorXd make_copy(std::vector<real> const& vec);

class Simulation {
public:

    Simulation(MeshGroup const& initial,real p_inf = 0.0, real epsilon = 1.0, real gamma = 1.0, real sigma = 0.0, real lambda = 1.0)
        :p_inf(p_inf), 
        epsilon(epsilon), 
        gamma(gamma), 
        sigma(sigma), 
        lambda(lambda),
        time(0.0),
        group(initial) {
            V_0 = group.volumes();

            inter.set_quadrature(quadrature_3);
            inter.set_quadrature(gauss_7);

#ifdef VERBOSE
            std::cout << "Simulation initialized with N = " << group.get_joined().verts.size() << " vertices and M = " << group.get_joined().trigs.size() << " triangles." << std::endl;
            std::cout << " epsilon: " << epsilon << std::endl;
            std::cout << " sigma:   " << sigma << std::endl;
            std::cout << " lambda:  " << lambda << std::endl;
#endif
        }

    virtual ~Simulation() {}

    virtual size_t phi_dim() const = 0;
    virtual size_t psi_dim() const = 0;

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const = 0;
    // if no other mesh is specified, the mesh describing the current state of the simulation is used:
    void assemble_matrices_prop(Eigen::MatrixXd& G,Eigen::MatrixXd& H) {
        assemble_matrices(G,H,group.get_joined());
    }

    Eigen::VectorXd solve_system(Eigen::MatrixXd const& G,Eigen::VectorXd const& H_phi) const;

    void compute_psi() {
        Eigen::MatrixXd G,H;
        assemble_matrices_prop(G,H);
        psi = solve_system(G,H*phi);
    }

    virtual void evolve_system(real dp, bool fixdt = false) = 0;

    void export_mesh(std::string fname);
    void export_mesh_values(std::string fname,std::vector<real> values);

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

    void set_V_0(std::vector<real> values) {
        V_0 = values;
    }

    real get_phi(size_t index) const {
        assert(index < phi_dim());
        return phi(index);
    }

    real get_psi(size_t index) const {
        assert(index < psi_dim());
        return psi(index);
    }

    real get_time() const {
        return time;
    }

    real get_volume() const {
        real vol(0.0);
        std::vector<real> vols(group.volumes());
        for(real elm : vols)
            vol += elm;
        return vol;
    }

    /*
    const std::vector<vec3>& get_vertices() {
        return mesh.verts;
    }*/

    // consider using psi.data() or similar instead of creating a new vector
    std::vector<real> get_psi() {
        return make_copy(psi);
    }

    std::vector<real> get_phi() {
        return make_copy(phi);
    }
    
protected:
    real get_dt(real dp,std::vector<vec3> const& gradients) const;
    real get_dt(real dp,std::vector<vec3> const& gradients, std::vector<real> const& kappa) const;
    real get_dt(real dp,std::vector<vec3> const& gradients, std::vector<real> const& kappa, real t) const;
    // (Material) time derivative of the potential as a function of the squared gradient.
    real potential_t(real grad_squared) const;
    real potential_t(real grad_squared,real kappa) const;
    real potential_t(real grad_squared, real volume, real kappa) const;
    real potential_t(real grad_squared, real volume, real kappa, real t) const;


    // The following constants are given in simulation units;
    // Eventual conversions happen outside this class

    // necessary constants (to run the simulation)
    real p_inf, epsilon, gamma, sigma, lambda;

    // initial condition informations:
    std::vector<real> V_0;

    // simulation time:
    real time;

    // Mesh describing the Bubble surface and an integrator object,
    // providing some functions to integrate over the mesh.
public:
    //Mesh mesh;
protected:
    Integrator inter;

    // Vectors that approximate the potential (phi) and its normal derivative (psi) on the surface
    Eigen::VectorXd psi;
    Eigen::VectorXd phi;

    MeshGroup group;
    

};

} // namespace Group

} // namespace Bem

#endif // SIMULATION_GROUP_HPP