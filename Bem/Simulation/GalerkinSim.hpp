#ifndef GALERKINSIM_HPP
#define GALERKINSIM_HPP

#include <iostream>
#include <vector>

#include "LinLinSim.hpp"

#include <Eigen/Dense>

namespace Bem {

class GalerkinSim : public LinLinSim {
public:

    GalerkinSim(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :LinLinSim(initial,p_inf,epsilon,sigma,gamma,pressurefield) {
            inter.set_quadrature(quadrature_7);
            inter.set_quadrature(gauss_7);
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    virtual ~GalerkinSim() {}

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;

};

} // namespace Bem

#endif // GALERKINSIM_HPP