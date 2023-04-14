#include "CurvatureTensor.hpp"

using namespace std;
using namespace Bem;
using Eigen::Matrix2d;
using Eigen::Vector;
using Eigen::Matrix;


void generate_curvature_tensor(Mesh const& mesh, vector<Bem::real>& kappa, vector<Bem::real>& gamma) {
    // initialize results with zero
    kappa = vector<Bem::real>(mesh.N(),0.0);
    gamma = vector<Bem::real>(mesh.N(),0.0);
    // also create a vector with weights
    vector<Bem::real> weights(mesh.N(),0.0);
    vector<vec3> vertex_normals = mesh.generate_vertex_normals_max();
    const vector<vec3>& vertices = mesh.get_vertices();
    vector<CoordSystem> vertex_systems(mesh.N());
    for(size_t i(0);i<mesh.N();++i) {
        vertex_systems[i] = CoordSystem(vertices[i],vertex_normals[i]);
        vertex_systems[i].swap_xz();
    }

    for(size_t i(0);i<mesh.M();++i) {
        Triplet t(mesh.get_triangle(i));
        vec3 e0 = vertices[t.b] - vertices[t.a];
        vec3 e1 = vertices[t.c] - vertices[t.b];
        vec3 e2 = vertices[t.a] - vertices[t.c];

        // build up coordinate system (u,v,n) on triangle
        vec3 u = e0;
        u.normalize();
        vec3 n = e0.vec(e1);
        n.normalize();
        vec3 v = n.vec(u);

        // compute the three components e,f,g of the curvature tensor
        // by using least square fitting of finite differences, see Rusinkiewicz_2004
        // the tensor is of the form | e  f | 
        //                           | f  g | and 0.5*(e+f) is the mean curvature

        Matrix<Bem::real,6,3> M {
            {e0.dot(u),e0.dot(v),0.0},
            {0.0,e0.dot(u),e0.dot(v)},
            {e1.dot(u),e1.dot(v),0.0},
            {0.0,e1.dot(u),e1.dot(v)},
            {e2.dot(u),e2.dot(v),0.0},
            {0.0,e2.dot(u),e2.dot(v)}
        };

        Vector<Bem::real,6> b {
            (vertex_normals[t.b]-vertex_normals[t.a]).dot(u),
            (vertex_normals[t.b]-vertex_normals[t.a]).dot(v),
            (vertex_normals[t.c]-vertex_normals[t.b]).dot(u),
            (vertex_normals[t.c]-vertex_normals[t.b]).dot(v),
            (vertex_normals[t.a]-vertex_normals[t.c]).dot(u),
            (vertex_normals[t.a]-vertex_normals[t.c]).dot(v)
        };

        Eigen::FullPivLU<Matrix<Bem::real,6,3>> solver;
        solver.compute(M);

        Vector<Bem::real,3> x = solver.solve(b);

        Bem::real mean_curvature = 0.5*(x(0)+x(2)); // Trace is base independent: we can take the trace of this matrix to obtain the mean curvature
        Bem::real gauss_curvature = x(0)*x(2) - x(1)*x(1); // determinant (base ind. too) is product of eigenvalues

        // Note we use here trivial weights, this could be improved on by taking for example
        // voronoi-like weightings (see Rusinkiewicz_2004)
        kappa[t.a] += mean_curvature;
        gamma[t.a] += gauss_curvature;
        weights[t.a] += 1.0;

        kappa[t.b] += mean_curvature;
        gamma[t.b] += gauss_curvature;
        weights[t.b] += 1.0;

        kappa[t.c] += mean_curvature;
        gamma[t.c] += gauss_curvature;
        weights[t.c] += 1.0;


        //cout << "mean_curvature: " << 0.5*(x(0)+x(2))<< endl;
        

    }

    for(size_t i(0);i<mesh.N();++i) {
        kappa[i] /= weights[i];
        gamma[i] /= weights[i];
    }
}