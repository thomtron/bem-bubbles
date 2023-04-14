#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "../basic/Bem.hpp"
#include "Mesh.hpp"
#include "FittingTool.hpp"

// This function computes the curvature tensor on each triangle and then 
// averages the values for mean and gaussian curvature on the mesh's vertices.
// These values are stored in kappa (mean) and gamma (gaussian) curvature.
// from these curvatures we can obtain the principal values of the curvature tensor.
// The used methods are taken from the paper: Rusinkiewicz_2004.pdf
void generate_curvature_tensor(Bem::Mesh const& mesh,std::vector<Bem::real>& kappa, std::vector<Bem::real>& gamma);

