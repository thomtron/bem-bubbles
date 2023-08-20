#include "Mesh.hpp"

#include <set>
#include <Eigen/Dense> // for curvature computation
#include "FittingTool.hpp" // for CoordSystem

using namespace std;


namespace Bem {

void Mesh::add(Mesh other,vec3 const& position) {
    size_t n_old = verts.size();
    for(vec3 const& vec : other.verts) {
        verts.push_back(vec+position);
    }
    for(Triplet t : other.trigs) {
        t.a += n_old;
        t.b += n_old;
        t.c += n_old;
        trigs.push_back(t);
    }
}
void Mesh::scale(real s) {
    for(vec3& elm : verts)
        elm *= s;
}

// not most efficient way, but we don't have to introduce matrices like that
void Mesh::rotate(vec3 v) {
    vec3 v_n = v;
    real theta = v.norm();
    v_n.normalize();
    for(vec3& pos : verts) {

        if(pos.vec(v_n).norm2() > 1e-20) {
            vec3 pos_para = pos.dot(v_n)*v_n;
            vec3 pos_orth = pos - pos_para;
            real a = pos_orth.norm();
            pos_orth *= (1.0/a);
            vec3 w = v_n.vec(pos_orth);

            pos = pos_para + a*(pos_orth*cos(theta)+w*sin(theta));
        }
    }
}

void Mesh::translate(vec3 v) {
    for(vec3& elm : verts)
        elm += v;
}

void Mesh::clear() {
    verts.clear();
    trigs.clear();
}

bool Mesh::check_validity() const {
    bool result = true;
    vector<bool> used(verts.size(),false);
    for(Triplet const& trig : trigs) {
        used[trig.a] = true;
        used[trig.b] = true;
        used[trig.c] = true;
    }
    for(bool elm : used) {
        if(not elm) result = false;
    }
    if(not result) {
        cout << "not all verts used!" << endl;
    }
    return result;
}

vector<vector<size_t>> generate_triangle_indices(Mesh const& mesh) {
    vector<vector<size_t>> triangle_indices(mesh.verts.size());
    size_t m(mesh.trigs.size());
    for(size_t i(0);i<m;++i){
        triangle_indices[mesh.trigs[i].a].push_back(i);
        triangle_indices[mesh.trigs[i].b].push_back(i);
        triangle_indices[mesh.trigs[i].c].push_back(i);
    }
    return triangle_indices;
}

vector<vector<size_t>> generate_2_ring(Mesh const& mesh, vector<vector<size_t>> const& neighbours) {
    vector<vector<size_t>> ring;

    size_t n(mesh.verts.size());
    for(size_t i(0);i<n;++i) {
        set<size_t> second_ring_indices;
        for(size_t k : neighbours[i]) {
            for(size_t l : neighbours[k]) {
                second_ring_indices.insert(l);
            }
        }

        ring.push_back(vector<size_t>(second_ring_indices.begin(),second_ring_indices.end()));
    }

    return ring;
}

vector<vector<size_t>> generate_2_ring(Mesh const& mesh) {
    return generate_2_ring(mesh,generate_neighbours(mesh));
}

vector<vector<size_t>> generate_neighbours(Mesh const& mesh) {
    size_t n(mesh.verts.size());
    vector<set<size_t>> neighbours(n);
    vector<vector<size_t>> result;
    for(Triplet t : mesh.trigs){
        neighbours[t.a].insert(t.b);
        neighbours[t.a].insert(t.c);

        neighbours[t.b].insert(t.c);
        neighbours[t.b].insert(t.a);

        neighbours[t.c].insert(t.a);
        neighbours[t.c].insert(t.b);
    }
    // copy the indices in the set into a vector array for convenience
    for(size_t i(0);i<n;++i) {
        result.push_back(vector<size_t>(neighbours[i].begin(),neighbours[i].end()));
    }

    return result;
}

vector<vec3> generate_triangle_normals  (Mesh const& mesh) {
    vector<vec3> triangle_normals;
    for(Triplet const& tri : mesh.trigs){
        vec3 normal((mesh.verts[tri.b]-mesh.verts[tri.a]).vec(mesh.verts[tri.c]-mesh.verts[tri.a]));
        normal.normalize();
        triangle_normals.push_back(normal);
    }
    return triangle_normals;
}

// approximation of vertex normals as proposed by Max_1999 (see Rusinkiewicz_2004). Is exact for vertices on a sphere
vector<vec3> generate_vertex_normals(Mesh const& mesh, vector<vector<size_t>> const& triangle_indices) {
    vector<vec3> vertex_normals;

    size_t n(mesh.verts.size());

    for(size_t i(0);i<n;++i) {
        vec3 normal;
        for(size_t j : triangle_indices[i]) {
            Triplet t(mesh.trigs[j]);
            t.cyclic_reorder(i);
            vec3 B(mesh.verts[t.b]-mesh.verts[t.a]);
            vec3 C(mesh.verts[t.c]-mesh.verts[t.a]);
            normal += B.vec(C)*(1.0/(B.norm2()*C.norm2()));
        }
        normal.normalize();
        vertex_normals.push_back(normal);
    }
    return vertex_normals;
}

vector<vec3> generate_vertex_normals(Mesh const& mesh) {
    return generate_vertex_normals(mesh,generate_triangle_indices(mesh));
}

// solid angle computation using spherical trigonometry (rule in the last line of this function, see Todhunter_1886)
real solid_angle_at_vertex(Mesh const& mesh, vector<vector<size_t>> const& triangle_indices, vector<vec3> const& triangle_normals, size_t i) {
    vector<Triplet> trigs;
    for(size_t j : triangle_indices[i]){
        Triplet trig(mesh.trigs[j]);
        trig.cyclic_reorder(i);
        trigs.push_back(trig);
    }
    vector<Triplet> sorted_trigs;
    vector<size_t> trig_indices;
    sorted_trigs.push_back(trigs[0]);
    trig_indices.push_back(triangle_indices[i][0]);
    size_t num(0);
    for(size_t j(1);j<trigs.size();++j){
        size_t new_ind = sorted_trigs[j-1].c;
        for(size_t k(0);k<trigs.size();++k){
            if(trigs[k].b == new_ind){
                sorted_trigs.push_back(trigs[k]);
                trig_indices.push_back(triangle_indices[i][k]);
                k = trigs.size(); // exit loop
                num++;
            }
        }
    }

    size_t n(sorted_trigs.size());
    real sol_ang(0.0);
    for(size_t j(0);j<n;++j){
        vec3 ray(mesh.verts[sorted_trigs[j].c]-mesh.verts[sorted_trigs[j].a]);
        ray.normalize();
        real angle = M_PI + asin((triangle_normals[trig_indices[j]].vec(triangle_normals[trig_indices[(j+1)%n]])).dot(ray));
        sol_ang += angle;
    }

    return sol_ang - (n-2)*M_PI;
}

real solid_angle_at_vertex(Mesh const& mesh, size_t i) {
    return solid_angle_at_vertex(mesh,generate_triangle_indices(mesh),generate_triangle_normals(mesh),i);
}

real volume(Mesh const& mesh) {
    real volume(0.0);

    for(Triplet t : mesh.trigs){
        vec3 a(mesh.verts[t.a]);
        vec3 b(mesh.verts[t.b]);
        vec3 c(mesh.verts[t.c]);

        volume += (a.vec(b)).dot(c)/6.0; // signed volume of tetraedron.
    }

    return volume;
}

vec3 centerofmass(Mesh const& mesh) {
    // to obtain the center of mass of the mesh, we'll compute the center of mass
    // of each tetraeder that is formed by a surface triangle and the origin. This
    // center of mass is just given by the arithmetic mean of the four edgepoints.
    // Then these vectors are summed up with their (signed) volume as weights.
    // Note that we assume a constant density in the mesh.
    vec3 center;
    real total_volume(0.0);

    for(Triplet t : mesh.trigs) {
        vec3 a(mesh.verts[t.a]);
        vec3 b(mesh.verts[t.b]);
        vec3 c(mesh.verts[t.c]);
        real volume = (a.vec(b)).dot(c)/6.0;

        center += 0.25*(a+b+c)*volume;
        total_volume += volume;
    }

    return center*(1.0/total_volume);
}

void to_centerofmass(Mesh& mesh) {
    vec3 center = centerofmass(mesh);
    for(vec3& elm : mesh.verts)
        elm -= center;
    return;
}



// This function computes the curvature tensor on each triangle and then 
// averages the values for mean and gaussian curvature on the mesh's vertices.
// These values are stored in kappa (mean) and gamma (gaussian) curvature.
// from these curvatures we can obtain the principal values of the curvature tensor.
// The used methods are taken from the paper: Rusinkiewicz_2004.pdf
void curvatures(Mesh const& mesh, vector<real>& kappa, vector<real>& gamma) {
    using Eigen::Matrix2d;
    using Eigen::Vector;
    using Eigen::Matrix;

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());
    vector<vec3> vertex_normals = generate_vertex_normals(mesh);
    const vector<vec3>& vertices = mesh.verts;
    vector<CoordSystem> vertex_systems(n);
    for(size_t i(0);i<n;++i) {
        vertex_systems[i] = CoordSystem(vertices[i],vertex_normals[i]);
        vertex_systems[i].swap_xz();
    }

    for(size_t i(0);i<m;++i) {
        Triplet t(mesh.trigs[i]);
        vec3 e0 = vertices[t.b] - vertices[t.a];
        vec3 e1 = vertices[t.c] - vertices[t.b];
        vec3 e2 = vertices[t.a] - vertices[t.c];

        // build up coordinate system (u,v,n) on triangle
        vec3 u = e0;
        u.normalize();
        vec3 normal = e0.vec(e1);
        normal.normalize();
        vec3 v = normal.vec(u);

        // compute the three components e,f,g of the curvature tensor
        // by using least square fitting of finite differences, see Rusinkiewicz_2004
        // the tensor is of the form | e  f | 
        //                           | f  g | and 0.5*(e+f) is the mean curvature

        Matrix<real,6,3> M {
            {e0.dot(u),e0.dot(v),0.0},
            {0.0,e0.dot(u),e0.dot(v)},
            {e1.dot(u),e1.dot(v),0.0},
            {0.0,e1.dot(u),e1.dot(v)},
            {e2.dot(u),e2.dot(v),0.0},
            {0.0,e2.dot(u),e2.dot(v)}
        };

        Vector<real,6> b {
            (vertex_normals[t.b]-vertex_normals[t.a]).dot(u),
            (vertex_normals[t.b]-vertex_normals[t.a]).dot(v),
            (vertex_normals[t.c]-vertex_normals[t.b]).dot(u),
            (vertex_normals[t.c]-vertex_normals[t.b]).dot(v),
            (vertex_normals[t.a]-vertex_normals[t.c]).dot(u),
            (vertex_normals[t.a]-vertex_normals[t.c]).dot(v)
        };

        Eigen::FullPivLU<Matrix<real,6,3>> solver;
        solver.compute(M);

        Vector<real,3> x = solver.solve(b);

        real mean_curvature = 0.5*(x(0)+x(2)); // Trace is base independent: we can take the trace of this matrix to obtain the mean curvature
        real gauss_curvature = x(0)*x(2) - x(1)*x(1); // determinant (base ind. too) is product of eigenvalues

        kappa.push_back(mean_curvature);
        gamma.push_back(gauss_curvature);
    }
    return;
}

vector<real> max_curvature(Mesh const& mesh) {
    // Note we use here trivial weights, this could be improved on by taking for example
    // voronoi-like weightings (see Rusinkiewicz_2004)
    vector<real> kappa, gamma;
    vector<real> weights(mesh.verts.size(),0.0);
    vector<real> max_curv(mesh.verts.size(),0.0);
    curvatures(mesh,kappa,gamma);
    for(size_t i(0);i<kappa.size();++i) {
        real mean_curvature = kappa[i];
        real gauss_curvature = gamma[i];
        real determinant = sqrt(mean_curvature*mean_curvature - gauss_curvature);
        real max_curvature = max(abs(mean_curvature+determinant),abs(mean_curvature-determinant));

        Triplet t(mesh.trigs[i]);

        max_curv[t.a] += max_curvature;
        weights[t.a] += 1.0;

        max_curv[t.b] += max_curvature;
        weights[t.b] += 1.0;

        max_curv[t.c] += max_curvature;
        weights[t.c] += 1.0;
    }

    for(size_t i(0);i<mesh.verts.size();++i) {
        max_curv[i] /= weights[i];
    }
    return max_curv;
}



// function for splitting/joining meshes

Mesh join_meshes(vector<Mesh> const& list) {
    Mesh joined;
    for(Mesh const& m : list) {
        joined.add(m);
    }
    return joined;
}

void update_connected(set<size_t>& connected,set<size_t>& front,vector<vector<size_t>> const& neighbours) {
    set<size_t> new_front;
    for(size_t elm : front) {
        for(size_t candidate : neighbours[elm]) {
            if(connected.insert(candidate).second) new_front.insert(candidate);
        }
    }
    front = new_front;
}

vector<Mesh> split_by_loose_parts(Mesh const& mesh) {
    vector<vector<size_t>> vert_perm;
    return split_by_loose_parts(mesh,vert_perm);
}

// split_by_loose_parts scans the mesh for parts that are not connected with each other.
// In our case we want to split up the mesh describing a group of bubbles into a vector
// of meshes describing each only one bubble.
vector<Mesh> split_by_loose_parts(Mesh const& mesh, vector<vector<size_t>>& vert_perm) {
    vector<Mesh> result;

    Mesh temp = mesh;
    vector<vector<size_t>> neighbours = generate_neighbours(temp);
    vector<vector<size_t>> trig_indices = generate_triangle_indices(temp);
    vector<bool> processed(temp.verts.size(),false);

    while(true) {

        size_t start = 0;
        while(start < processed.size() and processed[start] == true) {
            start++;
        }
        if(start == processed.size()) return result;

        set<size_t> connected_verts;
        set<size_t> front;
        connected_verts.insert(start);
        front.insert(start);

        while(front.size() > 0) update_connected(connected_verts,front,neighbours);

        set<size_t> connected_trigs;
        for(size_t vert : connected_verts)
            for(size_t trig : trig_indices[vert])
                connected_trigs.insert(trig);
        
        
        Mesh new_mesh;
        vector<size_t> inds;

        size_t i(0);
        for(size_t elm : connected_verts) {
            processed[elm] = true;
            // add the new vertices to the mesh 
            new_mesh.verts.push_back(temp.verts[elm]);
            inds.push_back(elm);

            // change the triplet's indices such that they are now pointing to the right part
            for(size_t trig : trig_indices[elm]) {
                Triplet& t(temp.trigs[trig]);
                t.cyclic_reorder(elm);
                t.a = i;
            }
            i++;
        }
        // add the new triangles to the mesh
        for(size_t elm : connected_trigs)
            new_mesh.trigs.push_back(temp.trigs[elm]);

        result.push_back(new_mesh);

        vert_perm.push_back(inds);
    
    }

    return result;
}


} // namespace Bem