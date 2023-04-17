#include "ResultTypes.hpp"

namespace Bem {

// these factors are used for the Galerkin simulation. These numbers are the integrals
// of the functhion 2pi (half of the solid angle), multiplied by the different types of 
// basis functions, over the unit triangle -> the sum of all parameters is thus pi, 
// since the unit triangle has area 1/2.

template<>
const real HomoPair<real>::identical_H_factor = M_PI;

template<>
const LinElm Pair<real,LinElm>::identical_H_factor = 
std::array<real,3>({
    M_PI/3.0,
    M_PI/3.0,
    M_PI/3.0
});

template<>
const LinLinElm HomoPair<LinLinElm>::identical_H_factor = 
std::array<real,9>({
    M_PI/6.0,
    M_PI/12.0,
    M_PI/12.0,
    M_PI/12.0,
    M_PI/6.0,
    M_PI/12.0,
    M_PI/12.0,
    M_PI/12.0,
    M_PI/6.0
});

} // namespace Bem