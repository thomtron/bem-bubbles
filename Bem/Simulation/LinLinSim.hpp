#ifndef LINLINSIM_HPP
#define LINLINSIM_HPP

#include <iostream>
#include <vector>
#include <limits>

#include "Simulation.hpp"

#include <Eigen/Dense>

#define LINEAR true
#define MIRROR_MESH true

namespace Bem {

PotVec compute_exterior_pot(CoordVec const& pos,Mesh M,PotVec phi,PotVec psi);

class LinLinSim : public Simulation {
public:

    LinLinSim(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :Simulation(initial,p_inf,epsilon,sigma,gamma,pressurefield),
        eps(1e-2),
        damping_factor(0.0),
        min_elm_size(0.0),
        max_elm_size(std::numeric_limits<real>::max()) {
            phi = Eigen::VectorXd::Zero(phi_dim());
            psi = Eigen::VectorXd::Zero(psi_dim());

            curvature_params = curvature_param();
            
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    virtual ~LinLinSim() {}

    virtual size_t phi_dim() const override {
        return mesh.verts.size();
    }
    virtual size_t psi_dim() const override {
        return mesh.verts.size();
    }
    
    virtual void evolve_system(real dp, bool fixdt = false) override;
    virtual void evolve_system_RK4(real dp, bool fixdt = false);

    virtual void remesh(real L);

    void set_damping_factor(real value) {
        damping_factor = value;
    }

    void set_minimum_element_size(real value) {
        min_elm_size = value;
    }

    void set_maximum_element_size(real value) {
        max_elm_size = value;
    }


    std::vector<real> kappa(Mesh const& m) const;

    virtual CoordVec position_t(Mesh const& m,PotVec& pot) const;
    PotVec   pot_t(Mesh const& m,CoordVec const& gradients, real t) const;
    PotVec   pot_t_multi(Mesh const& m,CoordVec const& gradients, real t) const;


    PotVec exterior_pot(CoordVec const& positions) const;

    CoordVec nopenetration(real eps, real dt,CoordVec const& x0,CoordVec& c) const;

protected:

    real eps;

    void test_negative() const;

    std::vector<vec3> generate_tangent_gradients(Mesh const& m, std::vector<real> const& pot) const;

    std::vector<real> curvature_param() const;

    real damping_factor;
    real min_elm_size, max_elm_size;
    PotVec curvature_params;
    
};

} // namespace Bem

#endif // LINLINSIM_HPP