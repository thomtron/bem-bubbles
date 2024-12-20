#ifndef COLOCSIM_HPP
#define COLOCSIM_HPP

#include <iostream>
#include <vector>

#include "LinLinSim.hpp"

#include <Eigen/Dense>

namespace Bem {

class ColocSimPin : public LinLinSim {
public:

    ColocSimPin(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :LinLinSim(initial,p_inf,epsilon,sigma,gamma,pressurefield),N_pin(0),index(0) {
            N_pin = set_x_boundary(mesh); // This is a problem! if one wants to initialize the simulation with phi and psi !!!

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

            // ignore the above, plans have shifted!
            
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    ColocSimPin(Mesh const& initial, PotVec phi_init, PotVec psi_init, real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0, real (*pressurefield)(vec3 x,real t) = &default_field)
        :LinLinSim(initial,p_inf,epsilon,sigma,gamma,pressurefield),N_pin(0),index(0) {
            N_pin = set_x_boundary(mesh,phi_init,psi_init);
            set_phi(phi_init);
            set_psi(psi_init);
            
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    virtual ~ColocSimPin() {}

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;

    virtual CoordVec position_t(Mesh const& m,PotVec& pot) const override;
    virtual void remesh(real L) override;
    //PotVec   pot_t(Mesh const& m,CoordVec const& gradients, real t) const;
    //PotVec   pot_t_multi(Mesh const& m,CoordVec const& gradients, real t) const;

    size_t N_pin;

    size_t index;
    
    void rearrange_boundary(Mesh& M,std::vector<bool> bound) const;
    size_t set_x_boundary(Mesh& M, PotVec& phi_t, PotVec& psi_t) const;
    size_t set_x_boundary(Mesh& M) const;

};

} // namespace Bem

#endif // COLOCSIM_HPP