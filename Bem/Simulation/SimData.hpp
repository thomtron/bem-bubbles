#ifndef SIMDATA_HPP
#define SIMDATA_HPP

#include <vector>
#include "../basic/Bem.hpp"

namespace Bem {

// object that holds the vertex positions and the potential values 
// for one Simulatoin. Will be used for the RK4 method.
// Note: possibility for optimization via lazy evaluation, but this
// functionalities probably aren't a bottleneck, so we let that to 
// a later time (possibly).
struct SimData {
    std::vector<vec3> pos;
    std::vector<real> phi;

    void operator*=(real scalar) {
        for(vec3& p : pos) { p*=scalar; }
        for(real& p : phi) { p*=scalar; }
    }

    void operator+=(SimData const& other) {
        size_t N(pos.size());
        for(size_t i(0);i<N;++i) { pos[i] += other.pos[i]; }
        N = phi.size();
        for(size_t i(0);i<N;++i) { phi[i] += other.phi[i]; }
    }
};

SimData operator*(SimData dat,real scalar);
SimData operator*(real scalar,SimData dat);

SimData operator+(SimData a,SimData const& b);

} // namespace Bem

#endif // SIMDATA_HPP