#ifndef CONLINGALERKINSIM_HPP
#define CONLINGALERKINSIM_HPP

#include <iostream>
#include <vector>

#include "Simulation.hpp"

#include <Eigen/Dense>

namespace Bem {

class ConLinGalerkinSim : public Simulation {
public:

    ConLinGalerkinSim(Mesh const& initial,real p_inf = 0.0, real epsilon = 1.0, real gamma = 1.0, real sigma = 0.0)
        :Simulation(initial,p_inf,epsilon,gamma,sigma) {
            phi = Eigen::VectorXd::Zero(phi_dim());
            psi = Eigen::VectorXd::Zero(psi_dim());
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and constant for psi." << std::endl;
#endif
        }

    virtual ~ConLinGalerkinSim() {}

    virtual size_t phi_dim() const override {
        return mesh.verts.size();
    }
    virtual size_t psi_dim() const override {
        return mesh.verts.size();
    }

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;
    virtual void evolve_system(real dp, bool fixdt = false) override {
        std::cout << dp << ',' << fixdt << std::endl; // prevent warnings
    }
    CoordVec position_t(Mesh const& m,PotVec const& pot) const;

private:
    CoordVec generate_tangent_gradients(Mesh const& m, PotVec const& pot) const;
};

} // namespace Bem

#endif // CONLINGALERKINSIM_HPP