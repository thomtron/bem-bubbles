#ifndef COLOCSIM_GROUP_HPP
#define COLOCSIM_GROUP_HPP

#include <iostream>
#include <vector>

#include "Simulation-group.hpp"
#include "SimData.hpp"

#include <Eigen/Dense>

namespace Bem {

namespace Group {

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

class ColocSim : public Simulation {
public:

    ColocSim(MeshGroup const& initial,real p_inf = 0.0, real epsilon = 1.0, real gamma = 1.0, real sigma = 0.0, real lambda = 1.0)
        :Simulation(initial,p_inf,epsilon,gamma,sigma,lambda) {
            phi = Eigen::VectorXd::Zero(phi_dim());
            psi = Eigen::VectorXd::Zero(psi_dim());
            damping_factor = 0.8;
            curvature_params = curvature_param();
#ifdef VERBOSE
            std::cout << "This simulation has linear elements for phi and constant for psi." << std::endl;
#endif
        }

    virtual ~ColocSim() {}

    virtual size_t phi_dim() const override {
        return group.num_verts();
    }
    virtual size_t psi_dim() const override {
        return group.num_verts();
    }

    virtual void assemble_matrices(Eigen::MatrixXd& G,Eigen::MatrixXd& H, Mesh const& m) const override;
    virtual void evolve_system(real dp, bool fixdt = false) override;
    void evolve_system_RK4(real dp, bool fixdt = false);

    SimData gradient(SimData const& X) const;

    std::vector<real> kappa(Mesh const& mesh) const;

    CoordVec position_t(Mesh const& mesh,PotVec const& pot) const;
    PotVec  pot_t(Mesh const& mesh,CoordVec const& gradients) const;

    std::vector<vec3> generate_gradients(std::vector<real>& kappa,Mesh const& m,Eigen::VectorXd const& phi,Eigen::VectorXd const& psi) const;
    std::vector<vec3> generate_gradients(std::vector<real>& kappa);
    std::vector<vec3> generate_tangent_gradients(Mesh const& m, std::vector<real> const& pot) const;

    std::vector<real> curvature_param();

    real damping_factor;
    PotVec curvature_params;
    void remesh(Mesh& mesh,std::vector<real>& pot,real L);
    void remesh(real L);

    void export_mesh_psi(std::string fname) const;
    void export_mesh(std::string fname) const;
};

} // namespace Group

} // namespace Bem

#endif // COLOCSIM_GROUP_HPP