#ifndef SQUARE_HPP
#define SQUARE_HPP

#include "Mesh.hpp"

namespace Bem {

class Square : public Mesh {
public:
    Square(unsigned int dimension) {
        create_unitsquare(dimension);
    }
    virtual ~Square() {}
private:
    void create_unitsquare(unsigned int dim);
};

} // namespace Bem

#endif // SQUARE_HPP