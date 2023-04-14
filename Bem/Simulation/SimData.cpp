#include "SimData.hpp"


namespace Bem {

SimData operator*(SimData dat,real scalar) {
    dat *= scalar;
    return dat;
}
SimData operator*(real scalar,SimData dat) {
    dat *= scalar;
    return dat;
}

SimData operator+(SimData a,SimData const& b) {
    a += b;
    return a;
}

} // namespace Bem