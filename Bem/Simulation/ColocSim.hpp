#ifndef COLOCSIM_HPP
#define COLOCSIM_HPP

#include <iostream>
#include <vector>

#include "LinLinSim.hpp"

#include <Eigen/Dense>

#define MIRROR_MESH true

namespace Bem {

class ColocSim : public LinLinSim {
public:

    ColocSim(Mesh const& initial,real p_inf = 1.0, real epsilon = 1.0, real sigma = 0.0, real gamma = 1.0,real (*pressurefield)(vec3 x,real t) = &default_field)
        :LinLinSim(initial,p_inf,epsilon,sigma,gamma,pressurefield) {
            
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and psi." << std::endl;
#endif
        }

    virtual ~ColocSim() {}

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;

};

} // namespace Bem

#endif // COLOCSIM_HPP