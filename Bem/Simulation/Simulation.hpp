#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <string>
#include <vector>

#include "../basic/Bem.hpp"
#include "../Mesh/Mesh.hpp"
#include "../Integration/Integrator.hpp"
#include "SimData.hpp"

#include <Eigen/Dense>

#include <cassert>

namespace Bem {

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

using CoordVec = std::vector<vec3>;
using PotVec = std::vector<real>;

template<typename T>
std::vector<T> operator*(real s,std::vector<T> vec) {
    for(T& elm : vec)
        elm *= s;
    return vec;
}
template<typename T>
std::vector<T> operator*(std::vector<T> vec, real s) {
    return s*vec;
}
template<typename T>
std::vector<T> operator+(std::vector<T> const& v1,std::vector<T> v2) {
    assert(v1.size() == v2.size());
    for(size_t i(0);i<v1.size();++i) {
        v2[i] += v1[i];
    }
    return v2;
}

std::vector<real> make_copy(Eigen::VectorXd const& vec);

Eigen::VectorXd make_copy(std::vector<real> const& vec);

real default_waveform(vec3 x,real t);

class Simulation {
public:

    Simulation(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*waveform)(vec3 x,real t) = &default_waveform)
        :p_inf(p_inf), 
        epsilon(epsilon),  
        sigma(sigma), 
        gamma(gamma),
        waveform(waveform),
        time(0.0),
        min_dt(-1.0),
        dp_balance(3.0),
        bicgstab(true),
        num_threads(100),
        mesh(initial) {
            V_0 = volume(initial);

            inter.set_quadrature(quadrature_7);
            inter.set_quadrature(gauss_7);

#ifdef VERBOSE
            std::cout << "Simulation initialized with N = " << mesh.verts.size() << " vertices and M = " << mesh.trigs.size() << " triangles." << std::endl;
            std::cout << " epsilon: " << epsilon << std::endl;
            std::cout << " sigma:   " << sigma << std::endl;
            std::cout << " gamma:  " << gamma << std::endl;
#endif
        }

    virtual ~Simulation() {}

    virtual size_t phi_dim() const = 0;
    virtual size_t psi_dim() const = 0;

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const = 0;
    // if no other mesh is specified, the mesh describing the current state of the simulation is used:
    void assemble_matrices_prop(Eigen::MatrixXd& G,Eigen::MatrixXd& H) const {
        assemble_matrices(G,H,mesh);
    }

    Eigen::VectorXd solve_system(Eigen::MatrixXd const& G,Eigen::VectorXd const& H_phi) const;

    void compute_psi() {
        Eigen::MatrixXd G,H;
        assemble_matrices_prop(G,H);
        psi = solve_system(G,H*phi);
    }

    virtual void evolve_system(real dp, bool fixdt = false) = 0;

    void export_mesh(std::string fname) const;
    void export_mesh_values(std::string fname,std::vector<real> values) const;

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

    real get_time() const {
        return time;
    }

    real get_volume() const {
        return volume(mesh);
    }

    
    const std::vector<vec3>& get_vertices() {
        return mesh.verts;
    }

    // consider using psi.data() or similar instead of creating a new vector
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
    real potential_t(real grad_squared, real volume, real kappa, vec3 pos, real t) const;


    // The following constants are given in simulation units;
    // Eventual conversions happen outside this class

    // necessary constants (to run the simulation)
    real p_inf, epsilon, sigma, gamma;

    // initial condition informations:
    real V_0;

    // wave informations:
    real (*waveform)(vec3 x,real t);

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