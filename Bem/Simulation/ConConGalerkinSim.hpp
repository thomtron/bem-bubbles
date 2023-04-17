#ifndef CONCONGALERKINSIM_HPP
#define CONCONGALERKINSIM_HPP

#include <iostream>

#include "Simulation.hpp"

#include <Eigen/Dense>

namespace Bem {

class ConConGalerkinSim : public Simulation {
public:

    ConConGalerkinSim(Mesh const& initial,real p_inf = 0.0, real epsilon = 1.0, real gamma = 1.0, real sigma = 0.0)
        :Simulation(initial,p_inf,epsilon,gamma,sigma) {
            phi = Eigen::VectorXd::Zero(phi_dim());
            psi = Eigen::VectorXd::Zero(psi_dim());
#ifdef VERBOSE
            std::cout << "This simulation has constant elements" << std::endl;
#endif
        }

    virtual size_t phi_dim() const override {
        return mesh.trigs.size();
    }
    virtual size_t psi_dim() const override {
        return mesh.trigs.size();
    }

    virtual ~ConConGalerkinSim() {}

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;
    virtual void evolve_system(real dp, bool fixdt = false) override {
        std::cout << dp << ',' << fixdt << std::endl; // prevent warnings
    }
    CoordVec position_t(Mesh const& m,PotVec const& pot) const;

private:
    CoordVec generate_tangent_gradients(Mesh const& m, PotVec const& pot) const;
};

} // namespace Bem

#endif // CONCONGALERKINSIM_HPP