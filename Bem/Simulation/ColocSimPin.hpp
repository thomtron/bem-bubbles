#ifndef COLOCSIM_HPP
#define COLOCSIM_HPP

#include <iostream>
#include <vector>

#include "LinLinSim.hpp"

#include <Eigen/Dense>

namespace Bem {

class ColocSimPin : public LinLinSim {
public:

    ColocSimPin(Mesh const& initial, size_t N_pin,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :LinLinSim(initial,p_inf,epsilon,sigma,gamma,pressurefield),N_pin(N_pin) {
            // The Simulation assumes that the last N_pin vertices of the Mesh initial are pinned, which means they do not change their position throughout
            // the simulation. Furthermore it is assumed that the vertices are pinned on an infinite flat plane; these vertices thus must lie in the same plane.
            // Only the first N-N_pin vertices act as colocation points such that the resulting linear system is again exactly determined (not overdetermined).
            // The matrix G has thus dimension N-N_pin x N-N_pin, the matrix H has dimensions N-N_pin x N. The remaining elements of the phi-vector are filled 
            // with the potential value of the wall (which is constant and equal to the potential at infinity). An additional term is added on the H side of the
            // equation, that accounts for the integration over the infinite plane (except for the part within the bubble). This integration however can be
            // conducted analytically by measuring the solid angle of this (perforated) plane, since phi is constant on the plane. For the G matrix these values
            // are not considered since psi is zero on the plane.

            // note that the potential at infinity also should come into the linear equation (see master thesis p. 15) - this is an error that was made when 
            // doing the numerical experiments with the oscillating bubbles! Have to investigate how this affects the results! However, the results look 
            // astonishingly similar to experimental results, which indicates that not a large deviation must be expected (or some better insights...)
            
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    virtual ~ColocSimPin() {}

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;

private:
    size_t N_pin;

};

} // namespace Bem

#endif // COLOCSIM_HPP