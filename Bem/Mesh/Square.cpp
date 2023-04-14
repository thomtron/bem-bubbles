#include "Square.hpp"

using namespace Bem;

void Square::create_unitsquare(unsigned int dim) {
    
    vertices.clear();
    triangles.clear();

    real dx = 1.0/(dim-1);

    for(unsigned int i(0);i<dim;++i){
        for(unsigned int j(0);j<dim;++j){
            vertices.push_back(vec3(i*dx,j*dx,0.0));
            if(i<dim-1 && j<dim-1) {
                triangles.push_back(Triplet(i*dim+j,(i+1)*dim+j,i*dim+(j+1)));
                triangles.push_back(Triplet((i+1)*dim+(j+1),i*dim+(j+1),(i+1)*dim+j));
            }
        }
    }

}