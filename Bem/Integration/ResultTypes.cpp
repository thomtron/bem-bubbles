#include "ResultTypes.hpp"

namespace Bem {

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