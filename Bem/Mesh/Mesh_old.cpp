#include "Mesh.hpp"
#include <string>
#include <vector>
#include <set>
#include <algorithm> // for sort()
#include <numeric>   // for iota()

using namespace std;

using namespace Bem;

// gives the solid angle of the surface from the viewpoint of the vertex i
Mesh::real Mesh::solid_angle_at_vertex(size_t i) const {
    vector<Triplet> trigs;
    vector<vector<size_t>> triangle_indices = generate_triangle_indices();
    vector<vec3> triangle_normals = generate_triangle_normals();
    for(size_t j : triangle_indices[i]){
        Triplet trig(triangles[j]);
        trig.cyclic_reorder(i);
        trigs.push_back(trig);

    }
    vector<Triplet> sorted_trigs;
    vector<size_t> trig_indices;
    sorted_trigs.push_back(trigs[0]);
    trig_indices.push_back(triangle_indices[i][0]);
    //cout << "size = " << trigs.size() << endl;
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
        vec3 ray(vertices[sorted_trigs[j].c]-vertices[sorted_trigs[j].a]);
        ray.normalize();
        real angle = M_PI + asin((triangle_normals[trig_indices[j]].vec(triangle_normals[trig_indices[(j+1)%n]])).dot(ray));
        sol_ang += angle;
    }

    return sol_ang - (n-2)*M_PI;
}


// volume() returns the inner volume of the mesh (the volume that the normals are pointing away from)
Mesh::real Mesh::volume() const {
    
    real volume(0.0);

    for(Triplet t : triangles){
        vec3 a(vertices[t.a]);
        vec3 b(vertices[t.b]);
        vec3 c(vertices[t.c]);

        volume += (a.vec(b)).dot(c)/6.0; // signed volume of tetraedron.
    }

    return volume;
}


vector<vector<size_t>> Mesh::generate_triangle_indices() const {
    vector<vector<size_t>> triangle_indices;
    size_t n(vertices.size());
    for(size_t i(0);i<n;++i){
        triangle_indices.push_back(vector<size_t>());
    }
    size_t m(triangles.size());
    for(size_t i(0);i<m;++i){
        triangle_indices[triangles[i].a].push_back(i);
        triangle_indices[triangles[i].b].push_back(i);
        triangle_indices[triangles[i].c].push_back(i);
    }
    return triangle_indices;
}

// containts also the vertex at the center
vector<vector<size_t>> Mesh::generate_2_ring() const {
    vector<vector<size_t>> ring;

    vector<vector<size_t>> neighbours = generate_neighbours();

    for(size_t i(0);i<N();++i) {
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


vector<vector<size_t>> Mesh::generate_neighbours() const {
    vector<set<size_t>> neighbours(N());
    vector<vector<size_t>> result;
    for(Triplet t : triangles){
        neighbours[t.a].insert(t.b);
        neighbours[t.a].insert(t.c);

        neighbours[t.b].insert(t.c);
        neighbours[t.b].insert(t.a);

        neighbours[t.c].insert(t.a);
        neighbours[t.c].insert(t.b);
    }
    // copy the indices in the set into a vector array for convenience
    for(size_t i(0);i<N();++i) {
        result.push_back(vector<size_t>(neighbours[i].begin(),neighbours[i].end()));
    }

    return result;
}

vector<vector<size_t>> Mesh::generate_edges() const {
    vector<vector<size_t>> result;

    vector<Tuplet> edges;
    vector<size_t> triangle_index;

    for(size_t i(0);i<M();++i) {
        Triplet trip(triangles[i]);

        edges.push_back(Tuplet(trip.a,trip.b));
        edges.push_back(Tuplet(trip.b,trip.c));
        edges.push_back(Tuplet(trip.c,trip.a));

        triangle_index.push_back(i);
        triangle_index.push_back(i);
        triangle_index.push_back(i);

    }

    // create an index list of the size of edges
    size_t K(edges.size());
    //cout << K << endl;
    vector<size_t> indices(K);
    iota(indices.begin(),indices.end(),0);
    // sort the indices (w.r.t. the alphabetical order of the edges -> remove duplicates)
    sort(indices.begin(),indices.end(),[&edges](size_t a,size_t b) { return edges[a] < edges[b]; });

    for(size_t i(0);i<K;i+=2) {
        size_t ind1 = indices[i];
        size_t ind2 = indices[i+1];
        vector<size_t> vals;
        vals.push_back(edges[ind1].get_a());
        vals.push_back(edges[ind1].get_b());
        vals.push_back(triangle_index[ind1]);
        vals.push_back(triangle_index[ind2]);
        result.push_back(vals);
        //cout << vals[0] << ' ' << vals[1] << ' ' << vals[2] << ' ' << vals[3] << endl;
    }

    return result;
}


#if 0
HalfedgeMesh Mesh::generate_halfedges() const {

    Halfedges result;

    // format halfedge: 
        //  - twin
        //  - next
        //  - vertex
        //  - edge
        //  - triangle
    
    // output arrays
    vector<size_t>& verts = result.verts;
    verts = vector<size_t>(N());
    vector<size_t>& faces = result.faces;
    vector<size_t>& edges = result.edges;
    vector<vector<size_t>>& halfs = result.halfs;

    // helper arrays
    vector<Tuplet> edge_verts;
    vector<size_t> triangle_index;



    for(size_t i(0);i<M();++i) {
        Triplet trip(triangles[i]);

        vector<size_t> halfedge(5);
        halfedge[1] = 3*i+1;
        halfedge[2] = trip.a;
        halfedge[4] = i;
        halfs.push_back(halfedge);

        halfedge[1] = 3*i+2;
        halfedge[2] = trip.b;
        halfs.push_back(halfedge);

        halfedge[1] = 3*i;
        halfedge[2] = trip.c;
        halfs.push_back(halfedge);


        faces.push_back(3*i);
        verts[trip.a] = 3*i;
        verts[trip.b] = 3*i+1;
        verts[trip.c] = 3*i+2;


        edge_verts.push_back(Tuplet(trip.a,trip.b));
        edge_verts.push_back(Tuplet(trip.b,trip.c));
        edge_verts.push_back(Tuplet(trip.c,trip.a));

        triangle_index.push_back(i);
        triangle_index.push_back(i);
        triangle_index.push_back(i);

    }

    // create an index list of the size of edges
    size_t K(edge_verts.size());
    //cout << K << endl;
    vector<size_t> indices(K);
    iota(indices.begin(),indices.end(),0);
    // sort the indices (w.r.t. the alphabetical order of the edge_verts -> remove duplicates)
    sort(indices.begin(),indices.end(),[&edge_verts](size_t a,size_t b) { return edge_verts[a] < edge_verts[b]; });


    for(size_t i(0);i<K;i+=2) {
        size_t ind1 = indices[i];
        size_t ind2 = indices[i+1];

        size_t trig1 = triangle_index[ind1];
        size_t trig2 = triangle_index[ind2];

        Triplet trip1 = triangles[trig1];
        Triplet trip2 = triangles[trig2];

        vector<size_t> shifts1;
        vector<size_t> shifts2;

        for(size_t s(0);s<3;++s){
            for(size_t t(0);t<3;++t){
                if(trip1[s] == trip2[t]) {
                    shifts1.push_back(s);
                    shifts2.push_back(t);
                }
            }
        }

        size_t shift1(shifts1[0]);
        size_t shift2(shifts2[0]);

        if((shifts1[1]+1)%3 == shifts1[0]) shift1 = shifts1[1];
        if((shifts2[1]+1)%3 == shifts2[0]) shift2 = shifts2[1];

        // push new halfedge index to edges array
        edges.push_back(3*trig1 + shift1);
        // write twin index
        halfs[3*trig1 + shift1][0] = 3*trig2 + shift2;
        halfs[3*trig2 + shift2][0] = 3*trig1 + shift1;
        // write edge index
        halfs[3*trig1 + shift1][3] = i;
        halfs[3*trig2 + shift2][3] = i;

        

    }



    /*
    for all edges:
        push back two half edges -> also push them back for 
        the faces 
        
        when generating edges we can easily generate halfedges
            - twin
            - next
            - vertex
            - trig
            - edge of course

            nooo twin is not that trivial! how do we know the facial neighbours of each triangle ?

            should have a function to loop over all direct-face-pairs of triangles! 

            we can generate the triangle indices for each vertex
            could basically do it like in the generate_edges function - not exactly!! ->
            make tuplets with triangle indices -> put em into a set

            then we have our triangle-pair-list!

            
        the other indices should be easily done too in the same for loop
        
    */

    return result;
}
#endif


HalfedgeMesh Mesh::generate_halfedges() const {

    HalfedgeMesh result;
    
    // output arrays
    vector<Halfedge*>& verts = result.verts;
    verts = vector<Halfedge*>(N(),nullptr);
    vector<Halfedge*>& faces = result.faces;
    vector<Halfedge*>& edges = result.edges;
    vector<Halfedge>& halfs = result.halfs;
    halfs = vector<Halfedge>(3*M());

    // helper arrays
    vector<Tuplet> edge_verts;
    vector<size_t> triangle_index;



    for(size_t i(0);i<M();++i) {
        Triplet trip(triangles[i]);

        halfs[3*i].next = &halfs[3*i+1];
        halfs[3*i].vertex = trip.a;
        halfs[3*i].triangle = i;

        halfs[3*i+1].next = &halfs[3*i+2];
        halfs[3*i+1].vertex = trip.b;
        halfs[3*i+1].triangle = i;

        halfs[3*i+2].next = &halfs[3*i];
        halfs[3*i+2].vertex = trip.c;
        halfs[3*i+2].triangle = i;


        faces.push_back(&halfs[3*i]);
        verts[trip.a] = &halfs[3*i];
        verts[trip.b] = &halfs[3*i+1];
        verts[trip.c] = &halfs[3*i+2];


        edge_verts.push_back(Tuplet(trip.a,trip.b));
        edge_verts.push_back(Tuplet(trip.b,trip.c));
        edge_verts.push_back(Tuplet(trip.c,trip.a));

        triangle_index.push_back(i);
        triangle_index.push_back(i);
        triangle_index.push_back(i);

    }

    // create an index list of the size of edges
    size_t K(edge_verts.size());
    //cout << K << endl;
    vector<size_t> indices(K);
    iota(indices.begin(),indices.end(),0);
    // sort the indices (w.r.t. the alphabetical order of the edge_verts -> remove duplicates)
    sort(indices.begin(),indices.end(),[&edge_verts](size_t a,size_t b) { return edge_verts[a] < edge_verts[b]; });


    for(size_t i(0);i<K/2;++i) {
        size_t ind1 = indices[2*i];
        size_t ind2 = indices[2*i+1];

        size_t trig1 = triangle_index[ind1];
        size_t trig2 = triangle_index[ind2];

        Triplet trip1 = triangles[trig1];
        Triplet trip2 = triangles[trig2];

        vector<size_t> shifts1;
        vector<size_t> shifts2;

        for(size_t s(0);s<3;++s){
            for(size_t t(0);t<3;++t){
                if(trip1[s] == trip2[t]) {
                    shifts1.push_back(s);
                    shifts2.push_back(t);
                }
            }
        }

        size_t shift1(shifts1[0]);
        size_t shift2(shifts2[0]);

        if((shifts1[1]+1)%3 == shifts1[0]) shift1 = shifts1[1];
        if((shifts2[1]+1)%3 == shifts2[0]) shift2 = shifts2[1];

        // push new halfedge index to edges array
        edges.push_back(&halfs[3*trig1 + shift1]);
        // write twin index
        halfs[3*trig1 + shift1].twin = &halfs[3*trig2 + shift2];
        halfs[3*trig2 + shift2].twin = &halfs[3*trig1 + shift1];
        // write edge index
        halfs[3*trig1 + shift1].edge = i;
        halfs[3*trig2 + shift2].edge = i;

        

    }

    return result;
}


vector<vector<size_t>> Mesh::generate_triangle_edge_indices(vector<vector<size_t>> const& edges) const {
    // initialize result with M elements
    vector<vector<size_t>> result(M());

    size_t K(edges.size());
    for(size_t i(0);i<K;++i){
        const vector<size_t>& edge(edges[i]);
        result[edge[2]].push_back(i);
        result[edge[3]].push_back(i);
    }

    return result;
}

void Mesh::split_edges(real L_max) {
    cout << "SPLIT-EDGES" << endl;
    vector<vector<size_t>> edges = generate_edges();
    vector<vector<size_t>> triangle_edges(generate_triangle_edge_indices(edges));
    vector<real> lengths;

    real L_max_2 = L_max*L_max;
    for(auto edge : edges) {
        lengths.push_back((vertices[edge[1]]-vertices[edge[0]]).norm2());
    }


    // we want to check each edge only once
    size_t J(edges.size());
    for(size_t i(0);i<J;++i) {

        vector<size_t> edge = edges[i];
        if(lengths[i] > L_max_2) {
            vertices.push_back(0.5*(vertices[edge[0]]+vertices[edge[1]]));
            size_t new_index(vertices.size()-1);

            size_t a(edge[2]);
            size_t b(edge[3]);
            
            Triplet A(triangles[a]);
            Triplet B(triangles[b]);
            A.cyclic_reorder(edge[0]);
            B.cyclic_reorder(edge[0]);

            if(not (A.b == edge[1])) {
                swap(a,b);
                swap(A,B);
            }

            size_t A_tip(0),B_tip(0);
            for(size_t k(0);k<3;++k) {
                if(A[k] != edge[0] and A[k] != edge[1]) A_tip = A[k];
                if(B[k] != edge[0] and B[k] != edge[1]) B_tip = B[k];
            }

            triangles[a] = Triplet(edge[0],new_index,A_tip);
            triangles[b] = Triplet(edge[0],B_tip,new_index);
            triangles.push_back( Triplet(new_index,edge[1],A_tip)); // num_trigs-2
            triangles.push_back( Triplet(new_index,B_tip,edge[1])); // num_trigs-1

            size_t num_trigs(triangles.size()); // for indexing

//------------------------------------------------------------------

            // getting to know which other edges belong to triangle A
            vector<size_t> tmp(triangle_edges[a]);
            Triplet edges_A(tmp[0],tmp[1],tmp[2]);
            edges_A.cyclic_reorder(i);

            // test whether edge edges_A.c is the one corresponding to edge[0]
            // -> ensure that edges_A.b is corresponding to edge[0]
            if(edges[edges_A.c][0] == edge[0] or edges[edges_A.c][1] == edge[0]) {
                swap(edges_A.c,edges_A.b);
            } 
            // for the edge .b we don't have to change anything.
            // for edge .c we have to adjust the triangle number to the newly generated one

            // check which of the two triangle indices corresponded to 'a' and set new coordinate
            if(edges[edges_A.c][3] == a) {
                edges[edges_A.c][3] = num_trigs - 2;
            } else {
                edges[edges_A.c][2] = num_trigs - 2;
            }

//------------------------------------------------------------------

            // getting to know which other edges belong to triangle B
            tmp = triangle_edges[b];
            Triplet edges_B(tmp[0],tmp[1],tmp[2]);
            edges_B.cyclic_reorder(i);

            // test whether edge edges_B.c is the one corresponding to edge[0]
            // -> ensure that edges_B.b is corresponding to edge[0]
            if(edges[edges_B.c][0] == edge[0] or edges[edges_B.c][1] == edge[0]) {
                swap(edges_B.c,edges_B.b);
            }
            // for the edge .b we don't have to change anything.
            // for edge .c we have to adjust the triangle number to the newly generated one

            // check which of the two triangle indices corresponded to 'b' and set new coordinate
            if(edges[edges_B.c][3] == b) {
                edges[edges_B.c][3] = num_trigs - 1;
            } else {
                edges[edges_B.c][2] = num_trigs - 1;
            }



            // Adding the newly created edges to the array "edges":
            // (we reset edges[i] (since the old one should be deleted)
            // and add the three newly created ones)
            // also compute the lengths of the new edges.

            vector<size_t> new_edge;
            Tuplet nedg(edge[0],new_index);
            new_edge.push_back(nedg.get_a());
            new_edge.push_back(nedg.get_b());
            new_edge.push_back(a);
            new_edge.push_back(b);

            edges[i] = new_edge; // i
            new_edge.clear();
            lengths[i] = (vertices[nedg.get_a()]-vertices[nedg.get_b()]).norm2();


            nedg = Tuplet(new_index,edge[1]);
            new_edge.push_back(nedg.get_a());
            new_edge.push_back(nedg.get_b());
            new_edge.push_back(num_trigs-2);
            new_edge.push_back(num_trigs-1);

            edges.push_back(new_edge); // num_edges-3
            new_edge.clear();
            lengths.push_back((vertices[nedg.get_a()]-vertices[nedg.get_b()]).norm2());


            nedg = Tuplet(A_tip,new_index);
            new_edge.push_back(nedg.get_a());
            new_edge.push_back(nedg.get_b());
            new_edge.push_back(num_trigs-2);
            new_edge.push_back(a);

            edges.push_back(new_edge); // num_edges-2
            new_edge.clear();
            lengths.push_back((vertices[nedg.get_a()]-vertices[nedg.get_b()]).norm2());


            nedg = Tuplet(B_tip,new_index);
            new_edge.push_back(nedg.get_a());
            new_edge.push_back(nedg.get_b());
            new_edge.push_back(b);
            new_edge.push_back(num_trigs-1);

            edges.push_back(new_edge); // num_edges-1
            new_edge.clear();
            lengths.push_back((vertices[nedg.get_a()]-vertices[nedg.get_b()]).norm2());



            // finally add the triangles to the triangle_indices list for the new edges:
            size_t num_edges(edges.size());

            triangle_edges[a].clear();
            triangle_edges[a].push_back(i);
            triangle_edges[a].push_back(num_edges-2);
            triangle_edges[a].push_back(edges_A.b);

            triangle_edges[b].clear();
            triangle_edges[b].push_back(i);
            triangle_edges[b].push_back(num_edges-1);
            triangle_edges[b].push_back(edges_B.b);

            // add two additional entries into the triangle_edges array
            triangle_edges.push_back(vector<size_t>());
            triangle_edges.push_back(vector<size_t>());
            
            triangle_edges[num_trigs-2].push_back(num_edges-3);
            triangle_edges[num_trigs-2].push_back(num_edges-2);
            triangle_edges[num_trigs-2].push_back(edges_A.c);

            triangle_edges[num_trigs-1].push_back(num_edges-3);
            triangle_edges[num_trigs-1].push_back(num_edges-1);
            triangle_edges[num_trigs-1].push_back(edges_B.c);

        }
    }
}


void remove_edge(vector<Halfedge*>& edges,size_t edge_ind) {
    if(edge_ind != edges.size()-1) {
        swap(edges[edge_ind],edges.back());
        Halfedge* u(edges[edge_ind]);
        u->edge = edge_ind;
        u->twin->edge = edge_ind;
    }
    edges.pop_back();
}

void remove_vertex(vector<Halfedge*>& verts,size_t vert_ind) {
    if(vert_ind != verts.size()-1) {
        swap(verts[vert_ind],verts.back());
        Halfedge* u(verts[vert_ind]);
        //cout << "vertex-pointer-II: " << u << endl;
        Halfedge* v(u);
        do {
            v->vertex = vert_ind;
            v = v->twin->next;
        } while(v != u);
    }
    verts.pop_back();
}

void remove_face(vector<Halfedge*>& faces,size_t face_ind) {
    if(face_ind != faces.size()-1) {
        swap(faces[face_ind],faces.back());
        Halfedge* u(faces[face_ind]);
        u->triangle = face_ind;
        u->next->triangle = face_ind;
        u->next->next->triangle = face_ind;
    }
    faces.pop_back();
}

void remove_halfedge(HalfedgeMesh& mesh,Halfedge* halfedge) {
    if(halfedge < &mesh.halfs.back()) {
        swap(*halfedge,mesh.halfs.back());
        halfedge->twin->twin           = halfedge;
        halfedge->next->next->next     = halfedge;
        mesh.verts[halfedge->vertex]   = halfedge;
        mesh.edges[halfedge->edge]     = halfedge;
        mesh.faces[halfedge->triangle] = halfedge;
        mesh.halfs.pop_back();
    }
}


void Mesh::collapse_edges(real L_min) {
    cout << "COLLAPSE-EDGES" << endl;

    size_t num_t_init = triangles.size();

    // shouldn't be done each time in practice!
    // could take Halfedges as argument for ex.
    HalfedgeMesh h(generate_halfedges());

    vector<real> lengths;
    real L_min_2 = L_min*L_min;
    
    for(auto halfedge : h.edges) {
        lengths.push_back((vertices[halfedge->next->vertex]-vertices[halfedge->vertex]).norm2());
    }

    // reorder the halfedge mesh and lengths such that shortest edges come first
    
    vector<size_t> indices(lengths.size());
    iota(indices.begin(),indices.end(),0);
    sort(indices.begin(),indices.end(),[&lengths](size_t a,size_t b) { return lengths[a] < lengths[b]; });

    vector<Halfedge*> edges_tmp = h.edges;
    vector<real> lengths_tmp = lengths;

    for(size_t i(0);i<lengths.size();++i) {
        size_t k = indices[i];
        h.edges[i] = edges_tmp[k];
        h.edges[i]->edge = i;
        h.edges[i]->twin->edge = i;
        lengths[i] = lengths_tmp[k];
        
    }

    for(size_t k(0);k<h.edges.size();++k) {

        if(lengths[k] < L_min_2) { // and not processed[k]) {

            // format halfedge: 
            // 0 - twin
            // 1 - next
            // 2 - vertex
            // 3 - edge
            // 4 - triangle

            Halfedge* halfedge_A = h.edges[k];
            Halfedge* halfedge_B = halfedge_A->twin; // twin

            size_t vert_0 = halfedge_A->vertex;
            size_t vert_1 = halfedge_B->vertex; // next

            size_t vert_A = halfedge_A->next->next->vertex;
            size_t vert_B = halfedge_B->next->next->vertex;

            Halfedge* half_A = h.verts[vert_A];
            Halfedge* current = half_A->twin->next;
            size_t valence_A = 1;
            while( current != half_A) {
                valence_A++;
                current = current->twin->next;
            }

            Halfedge* half_B = h.verts[vert_B];
            current = half_B->twin->next;
            size_t valence_B = 1;
            while( current != half_B) {
                valence_B++;
                current = current->twin->next;
            }

            // ensure that we wont produce "triangle appendices"
            if(valence_A > 3 and valence_B > 3) {
                // test whether all triangles are still
                // facing (approximately) in the right direction

                vec3 new_pos(0.5*(vertices[vert_0] + vertices[vert_1]));
                vec3 old_pos(vertices[vert_1]);

                bool valid = true;

                double limit_value = 0.8;

                Halfedge* u(halfedge_A->next->twin);
                do {
                    if(u == halfedge_B->next->next) {
                        u = halfedge_B->next->twin;
                        old_pos = vertices[vert_0];
                    }
                    vec3 a = vertices[u->vertex] - vertices[u->next->next->vertex];
                    vec3 b = vertices[u->next->next->vertex] - new_pos;
                    vec3 c = vertices[u->next->next->vertex] - old_pos;

                    c = c.vec(a);
                    c.normalize();
                    b = b.vec(a);
                    b.normalize();

                    if(b.dot(c)<limit_value) valid = false;
                    u = u->next->twin;

                } while(u != halfedge_A->next->next);


                if(valid) {

                    // remove the things

                    // remove the three edges
                    remove_edge(h.edges,halfedge_A->next->edge);
                    remove_edge(h.edges,halfedge_B->next->next->edge);
                    remove_edge(h.edges,halfedge_A->edge);

                    swap(lengths[halfedge_A->next->edge],lengths.back());
                    lengths.pop_back();
                    swap(lengths[halfedge_B->next->next->edge],lengths.back());
                    lengths.pop_back();
                    swap(lengths[halfedge_A->edge],lengths.back());
                    lengths.pop_back();


                    // remove the two faces
                    remove_face(h.faces,halfedge_A->triangle);
                    remove_face(h.faces,halfedge_B->triangle);


                    // make sure that vert_0 has an halfedge index that wont be removed
                    // the same for vert_A and vert_B
                    h.verts[vert_0] = halfedge_A->next->next->twin;
                    h.verts[vert_A] = halfedge_A->next->twin;
                    h.verts[vert_B] = halfedge_B->next->twin;

                    remove_vertex(h.verts,vert_1);

                    // adjust the real vertices of the mesh
                    vertices[vert_0] = new_pos;
                    swap(vertices[vert_1],vertices.back());
                    vertices.pop_back();

                    // in case that we just swapped vert_0, adjust it here (since we use it later on)
                    if(vert_0 == h.verts.size()) vert_0 = vert_1;

                    // move pointers from vert_1 to vert_0
                    current = halfedge_A->next;
                    while( current != halfedge_B) {
                        current->vertex = vert_0;
                        current = current->twin->next;
                    }

                    // make sure that the edges of triangle A and B that WONT be removed, 
                    // will still have a valid halfedge pointer
                    h.edges[halfedge_A->next->next->edge] = halfedge_A->next->next->twin;
                    h.edges[halfedge_B->next->edge] = halfedge_B->next->twin;

                
                    // first handle twins (before we destroy structure)
                    halfedge_A->next->next->twin->twin = halfedge_A->next->twin;
                    halfedge_A->next->twin->twin = halfedge_A->next->next->twin;
                    halfedge_A->next->twin->edge = halfedge_A->next->next->edge;

                    halfedge_B->next->next->twin->twin = halfedge_B->next->twin;
                    halfedge_B->next->twin->twin = halfedge_B->next->next->twin;
                    halfedge_B->next->next->twin->edge = halfedge_B->next->edge;

                    // finally remove the halfedges

                    halfedge_A->vertex = h.verts.size();
                    halfedge_A->next->vertex = h.verts.size();
                    halfedge_A->next->next->vertex = h.verts.size();
                    halfedge_B->vertex = h.verts.size();
                    halfedge_B->next->vertex = h.verts.size();
                    halfedge_B->next->next->vertex = h.verts.size();

                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_A->next->next);
                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_A->next);
                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_A);
                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_B->next->next);
                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_B->next);
                    while(h.halfs.back().vertex == h.verts.size()) h.halfs.pop_back();
                    remove_halfedge(h,halfedge_B);


                    // adjust lengths and valences

                    Halfedge* half_0 = h.verts[vert_0];
                    current = half_0;
                    do {
                        lengths[current->edge] = (vertices[current->next->vertex]-vertices[current->vertex]).norm2();
                        current = current->twin->next;
                    } while( current != half_0);

                    
                    k = 0;
                }
            }
        }
    }

    triangles.clear();
    for(Halfedge* half : h.faces) {
        triangles.push_back(Triplet(half->vertex,half->next->vertex,half->next->next->vertex));
    }

    cout << endl << "simplified mesh from " << num_t_init << " to " << triangles.size() << " triangles." << endl;

}

/* could be useful above...
size_t get_valence(Halfedge* start) {
    Halfedge* v = start;
    size_t valence(0);
    do {
        valence++;
        v = v->twin->next;
    } while(v != start);
    return valence;
}*/

int valence_cost(int valence) {
    int diff(valence - 6); // 6 is the optimal valence number
    return diff*diff;
}

vec3 cross_vec(vec3 const& a,vec3 const& b,vec3 const& c) {
    return (b-a).vec((c-b));
}

void Mesh::flip_edges(int state) {
    cout << "FLIP-EDGES" << endl;

    HalfedgeMesh mesh(generate_halfedges());
    bool v_cost = true;
    bool s_cost = true;
    if(state == 1) s_cost = false;
    if(state == 2) v_cost = false;

    vector<size_t> valences(N(),0);

    for(Triplet tri : triangles) {
        valences[tri.a]++;
        valences[tri.b]++;
        valences[tri.c]++;
    }
    

    size_t num_flipped(0);

    for(Halfedge* halfedge_A : mesh.edges){

        Halfedge* halfedge_B = halfedge_A->twin;

        size_t vert_0 = halfedge_A->vertex;
        size_t vert_1 = halfedge_B->vertex;

        size_t vert_A = halfedge_A->next->next->vertex;
        size_t vert_B = halfedge_B->next->next->vertex;
        


        if(valences[vert_0] > 3 and valences[vert_1] > 3) {

            int cost = 0;
            int cost_flip = 0;
            if(v_cost) {

                cost += valence_cost(valences[vert_0]);
                cost += valence_cost(valences[vert_1]);
                cost += valence_cost(valences[vert_A]);
                cost += valence_cost(valences[vert_B]);

                cost_flip += valence_cost(valences[vert_0]-1); // it is safe to subtract one here
                cost_flip += valence_cost(valences[vert_1]-1);
                cost_flip += valence_cost(valences[vert_A]+1);
                cost_flip += valence_cost(valences[vert_B]+1);

            }


            // introduce surface cost!!
            real surf_cost = 0.0;
            real surf_cost_flip = 0.0;

            if(s_cost) {
                surf_cost += cross_vec(vertices[vert_0],vertices[vert_1],vertices[vert_A]).norm2();
                surf_cost += cross_vec(vertices[vert_0],vertices[vert_B],vertices[vert_1]).norm2();

                surf_cost_flip += cross_vec(vertices[vert_0],vertices[vert_B],vertices[vert_A]).norm2();
                surf_cost_flip += cross_vec(vertices[vert_B],vertices[vert_1],vertices[vert_A]).norm2();

                
            }

            bool flip(true);
            if(v_cost) flip = (flip and cost_flip <= cost); // and surf_cost_flip/surf_cost < 1.2);
            if(s_cost) flip = (flip and surf_cost_flip < surf_cost);



            // maybe test also for delaunay condition ?? 

            vec3 nA = cross_vec(vertices[vert_0],vertices[vert_B],vertices[vert_A]);
            nA.normalize();
            vec3 nB = cross_vec(vertices[vert_B],vertices[vert_1],vertices[vert_A]);
            nB.normalize();
            real dotprod = nA.dot(nB);

            if(flip and dotprod > 0.95) { 
                valences[vert_0]--;
                valences[vert_1]--;
                valences[vert_A]++;
                valences[vert_B]++;

                // flip the edge!

                size_t edge = halfedge_A->edge;

                halfedge_A->edge = halfedge_B->next->edge;
                halfedge_A->twin = halfedge_B->next->twin;
                halfedge_B->next->twin->twin = halfedge_A;

                halfedge_B->edge = halfedge_A->next->edge;
                halfedge_B->twin = halfedge_A->next->twin;
                halfedge_A->next->twin->twin = halfedge_B;

                halfedge_A->next->vertex = vert_B;
                halfedge_A->next->twin = halfedge_B->next;
                halfedge_A->next->edge = edge;

                halfedge_B->next->vertex = vert_A;
                halfedge_B->next->twin = halfedge_A->next;
                halfedge_B->next->edge = edge;

                mesh.edges[edge] = halfedge_A->next;
                mesh.edges[halfedge_A->edge] = halfedge_A;
                mesh.edges[halfedge_B->edge] = halfedge_B;

                num_flipped++;
                
            }

        } 
        
    }

    triangles.clear();
    for(Halfedge* half : mesh.faces) {
        triangles.push_back(Triplet(half->vertex,half->next->vertex,half->next->next->vertex));
    }
    cout << "flipped " << num_flipped << " edges." << endl;
}

// need to reproject the vertices on to the original surface
void Mesh::relax_vertices() {
    cout << "RELAX-VERTICES" << endl;
    
    vector<vector<size_t>> neighbours = generate_neighbours();
    vector<vec3> new_vertices(N());

    for(size_t i(0);i<N();++i) {
        vec3 new_pos;
        real num(0.0);
        for(size_t j : neighbours[i]) {
            new_pos += vertices[j];
            num++;
        }
        new_vertices[i] = (1.0/num)*new_pos;
    }

    vertices = new_vertices;

}

#include <omp.h>


void Mesh::trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index) {     
    
    real s_min = -1.0;
    //cout << "Start:" << endl;

    for(size_t j(0);j<mesh.M();++j) {
        Triplet t(mesh.triangles[j]);

        vec3 a(mesh.vertices[t.b]-mesh.vertices[t.a]);
        vec3 b(mesh.vertices[t.c]-mesh.vertices[t.b]);
        vec3 c(mesh.vertices[t.a]-mesh.vertices[t.c]);

        vec3 n(a.vec(b));

        real s = n.dot(mesh.vertices[t.a]-pos)/n.dot(dir);

        vec3 x = pos + s*dir;

        // find the position with minimal s - parameter
        if(abs(s) < s_min or s_min < 0.0) { 
            if(     (mesh.vertices[t.a]-x).vec(a).dot(n) >= 0
                and (mesh.vertices[t.b]-x).vec(b).dot(n) >= 0
                and (mesh.vertices[t.c]-x).vec(c).dot(n) >= 0 ) {
                    //cout << "position found for " << i << "!" << endl;
                    s_min = abs(s);
                    result = x;
                    trig_index = j;
            }
        }
    }
}

void Mesh::project(Mesh const& other) {
    // normalized vertex normals
    vector<vec3> normals = generate_vertex_normals_max();

    vector<vec3> new_vertices(N());

    size_t n(N());

    omp_set_num_threads(100);

    #pragma omp parallel
    {

    Mesh local_this,local_other;
    vector<vec3> local_normals;

    #pragma omp critical 
    {
        local_this = *this;
        local_other = other;
        local_normals = normals;
    }


    #pragma omp for
    for(size_t i = 0;i<n;++i) {

        vec3 projected_pos;

        size_t index(0);

        trace_mesh(local_other,local_this.vertices[i],local_normals[i],projected_pos,index);

        #pragma omp critical
        new_vertices[i] = projected_pos;
    }

    }

    vertices.clear();
    vertices = new_vertices;
}


#include "FittingTool.hpp"

vector<Bem::real> Mesh::project_and_interpolate(Mesh const& other,vector<Bem::real> const& f) {
    cout << "PROJECT-AND-INTERPOLATE" << endl;

    assert(f.size() == other.N());
    vector<Bem::real> result(N());

    // normalized vertex normals
    vector<vec3> normals = generate_vertex_normals_max();

    vector<vec3> new_vertices(N());

    vector<vec3> other_normals = other.generate_vertex_normals_max();
    vector<vector<size_t>> verts = other.generate_neighbours(); //generate_2_ring();
    vector<FittingTool> fits(other.N());
    for(size_t i(0);i<other.N();++i) {
        vector<vec3> positions;
        positions.push_back(other.vertices[i]);
        for(size_t j : verts[i])
            positions.push_back(other.vertices[j]);
        fits[i].compute_quadratic_fit(other_normals[i],other.vertices[i],positions);
    }


    omp_set_num_threads(100);

    #pragma omp parallel
    {

    Mesh local_this,local_other;
    vector<vec3> local_normals;

    #pragma omp critical 
    {
        local_this = *this;
        local_other = other;
        local_normals = normals;
    }

    size_t n(local_this.N());

    #pragma omp for
    for(size_t i = 0;i<n;++i) {

        vec3 projected_pos;

        size_t index(0);

        trace_mesh(local_other,local_this.vertices[i],local_normals[i],projected_pos,index);

        Triplet t = local_other.get_triangle(index);

        // copy coord system of the three fits at the three edges.
        // using transform, and then get_position, obtain the projection on the 
        // current fit. then obtain distance of point to edges -> use linear weights
        // for averaging the positions obtained through the three fits.

        FittingTool A,B,C;
        #pragma omp critical 
        {
            A = fits[t.a];
            B = fits[t.b];
            C = fits[t.c];
        }

        CoordSystem Sa,Sb,Sc;
        Sa = A.copy_coord_system();
        Sb = B.copy_coord_system();
        Sc = C.copy_coord_system();

        vec3 Pa,Pb,Pc;
        Pa = Sa.transform(projected_pos);
        Pb = Sb.transform(projected_pos);
        Pc = Sc.transform(projected_pos);

        Pa = A.get_position(Pa.x,Pa.y);
        Pb = B.get_position(Pb.x,Pb.y);
        Pc = C.get_position(Pc.x,Pc.y);

        // now determine distances

        // build local 2d coord system
        vec3 u(local_other.vertices[t.b]-local_other.vertices[t.a]);
        u.normalize();
        vec3 v(local_other.vertices[t.c]-local_other.vertices[t.b]);
        v = v - v.dot(u)*u;
        v.normalize();

        vec3 a(local_other.vertices[t.b]-local_other.vertices[t.a]);
        vec3 b(local_other.vertices[t.c]-local_other.vertices[t.b]);
        vec3 x(projected_pos - local_other.vertices[t.a]);

        Bem::real a0,a1,b0,b1,x0,x1,q,r;

        a0 = a.dot(u);
        a1 = a.dot(v);
        b0 = b.dot(u);
        b1 = b.dot(v);
        x0 = x.dot(u);
        x1 = x.dot(v);

        q = (x1*b0 - b1*x0)/(a1*b0-a0*b1);
        r = (a1*x0 - a0*x1)/(a1*b0-a0*b1);

        projected_pos = (1.0 - q)*Pa + (q - r)*Pb + r*Pc;

        #pragma omp critical 
        {
            new_vertices[i] = projected_pos;
            result[i] = (1.0 - q)*f[t.a] + (q - r)*f[t.b] + r*f[t.c];
        }
    }

    }

    vertices.clear();
    vertices = new_vertices;

    return result;
}




vector<vec3> Mesh::generate_triangle_normals() const {
    vector<vec3> triangle_normals;
    for(Triplet const& tri : triangles){
        // calculate here the outward facing normal
        vec3 normal((vertices[tri.b]-vertices[tri.a]).vec(vertices[tri.c]-vertices[tri.a]));
        normal.normalize();
        triangle_normals.push_back(normal);
    }
    return triangle_normals;
}


// approximation of vertex normals as proposed by Max_1999. Is exact for vertices on a sphere
vector<vec3> Mesh::generate_vertex_normals_max() const {
    vector<vec3> vertex_normals;
    vector<vector<size_t>> triangle_indices = generate_triangle_indices();

    for(size_t i(0);i<N();++i) {
        vec3 normal;
        for(size_t j : triangle_indices[i]) {
            Triplet t(triangles[j]);
            t.cyclic_reorder(i);
            vec3 B(vertices[t.b]-vertices[t.a]);
            vec3 C(vertices[t.c]-vertices[t.a]);
            normal += B.vec(C)*(1.0/(B.norm2()*C.norm2()));
        }
        normal.normalize();
        vertex_normals.push_back(normal);
    }
    return vertex_normals;
}


// including up to second nearest neighbouring triangles
// note this still isn't according perfectly to ning's definition
vector<vec3> Mesh::generate_vertex_normals() const {

    vector<vec3> vertex_normals;
    vector<vector<size_t>> triangle_indices = generate_triangle_indices();
    vector<vec3> triangle_normals = generate_triangle_normals();

    for(size_t i(0);i<triangle_indices.size();++i){
        
        set<size_t> vert_indices;
        for(size_t ind : triangle_indices[i]){           
            Triplet tri(triangles[ind]);
            vert_indices.insert(tri.a);
            vert_indices.insert(tri.b);
            vert_indices.insert(tri.c);
        }
        set<size_t> tria_indices;
        for(size_t vert : vert_indices){
            for(size_t ind : triangle_indices[vert]){
                tria_indices.insert(ind);
            }
        }

        vector<real> dists;
        real mean_dist(0.0);
        for(size_t ind : tria_indices){
            Triplet t(triangles[ind]);
            real dist = (vertices[i] - (vertices[t.a] + vertices[t.b] + vertices[t.c]) * (1.0/3.0) ).norm();
            mean_dist += dist;
            dists.push_back(dist);
        }
        mean_dist /= real(dists.size());

        vec3 normal;
        size_t k(0);
        for(size_t ind : tria_indices){
            normal += triangle_normals[ind] * exp(-dists[k]/(2.0*mean_dist));
            k++;
        }

        normal.normalize();
        vertex_normals.push_back(normal);
    }

    return vertex_normals;
}

void Mesh::add_mesh(Mesh const& other,vec3 position) {
    size_t n_old = vertices.size();
    for(vec3 const& vec : other.vertices) {
        vertices.push_back(vec+position);
    }
    for(Triplet t : other.triangles) {
        t.a += n_old;
        t.b += n_old;
        t.c += n_old;
        triangles.push_back(t);
    }
}

vector<Mesh> Mesh::generate_mesh_partition(size_t num) const {
    vector<Mesh> partition;
    size_t num_elms = M()/num;

    if(num_elms == 0) {
        partition.push_back(*this);
        return partition;
    }

    size_t remainder = M()%num;

    size_t offset = 0;
    for(size_t i(0);i<remainder;++i) {
        Mesh m;
        m.vertices = vertices;
        for(size_t j(0);j<num_elms+1;++j) {
            m.add_triangle(triangles[j+offset]);
        }
        partition.push_back(m);
        offset += num_elms+1;
    }
    for(size_t i(0);i<num-remainder;++i) {
        Mesh m;
        m.vertices = vertices;
        for(size_t j(0);j<num_elms;++j) {
            m.add_triangle(triangles[j+offset]);
        }
        partition.push_back(m);
        offset += num_elms;
    }

    return partition;
}





#include <iostream>
#include <fstream> // for ofstream
#include <algorithm> // for clamp

void Mesh::export_obj(string filename) const
{
    ofstream output(filename);

    output << "# output from Mesh class" << endl;
    
    for(vec3 const& v : vertices){
        output << "v " << v.x << ' ' << v.y << ' ' << v.z << endl;
    }

    // note in .obj the indices start with 1 - we thus have to add 1 to each index.
    for(Triplet const& t : triangles){
        output << "f " << t.a+1 << ' ' << t.b+1 << ' ' << t.c+1 << endl;
    }

    output.close();
}


// these are two small structs needed for export_ply, since we want to export in binary form and 
// export floats and ints always in 32 bits (eventually contrary to the representations used in computations)
struct float_vec {
    float x,y,z;
    float_vec(float x,float y,float z)
        :x(x),y(y),z(z) {}
};

struct uint_vec {
    unsigned int i,j,k;
    uint_vec(unsigned int i,unsigned int j,unsigned int k)
        :i(i),j(j),k(k) {}
};

template<typename T>
struct Vertex {
    T x,y,z,w;
    Vertex(T x,T y,T z,T w)
        :x(x),y(y),z(z),w(w) {}
};

void Mesh::export_ply(string filename) const {
    ofstream output(filename);

    size_t n(N());
    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property float x" << endl; 
    output << "property float y" << endl; 
    output << "property float z" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        float_vec vec(vertices[i].x,vertices[i].y,vertices[i].z);
        output.write((char*) (&vec),sizeof(float_vec));
    }
    for(Triplet const& t : triangles){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void Mesh::export_ply(string filename,vector<real> values,real min,real max) const {
    ofstream output(filename);

    size_t n(N());
    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property float x" << endl; 
    output << "property float y" << endl; 
    output << "property float z" << endl;   
    output << "property uchar red"   << endl; 
    output << "property uchar green" << endl; 
    output << "property uchar blue"  << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        float_vec vec(vertices[i].x,vertices[i].y,vertices[i].z);
        output.write((char*) (&vec),sizeof(float_vec));
        char red(clamp((values[i]-min)/(max-min)*255,0.,255.));
        char green(0);
        char blue(clamp((max-values[i])/(max-min)*255,0.,255.));
        output.write(&red,sizeof(char));
        output.write(&green,sizeof(char));
        output.write(&blue,sizeof(char));
    }
    for(Triplet const& t : triangles){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}


void Mesh::export_ply_float(string filename,vector<real> values) const {
    ofstream output(filename);

    size_t n(N());
    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property float x" << endl; 
    output << "property float y" << endl; 
    output << "property float z" << endl;
    output << "property float w" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        Vertex<float> vec(vertices[i].x,vertices[i].y,vertices[i].z,values[i]);
        output.write((char*) (&vec),sizeof(Vertex<float>));
    }
    for(Triplet const& t : triangles){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void Mesh::export_ply_double(string filename,vector<real> values) const {
    ofstream output(filename);

    size_t n(N());
    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property double x" << endl; 
    output << "property double y" << endl; 
    output << "property double z" << endl;
    output << "property double w" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        Vertex<double> vec(vertices[i].x,vertices[i].y,vertices[i].z,values[i]);
        output.write((char*) (&vec),sizeof(Vertex<double>));
    }
    for(Triplet const& t : triangles){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void Mesh::export_ply_float_separat(string filename,vector<real> values) const {
    ofstream output(filename);

    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << 3*m << endl;
    output << "property float x" << endl; 
    output << "property float y" << endl; 
    output << "property float z" << endl;
    output << "property float w" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<m;++i) {
        for(size_t j(0);j<3;++j){
            size_t k = triangles[i][j];
            Vertex<float> vec(vertices[k].x,vertices[k].y,vertices[k].z,values[3*i+j]);
            output.write((char*) (&vec),sizeof(Vertex<float>));
        }
    }
    for(size_t i(0);i<m;++i){
        uint_vec vec(3*i,3*i+1,3*i+2);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void Mesh::export_ply_double_separat(string filename,vector<real> values) const {
    ofstream output(filename);

    size_t m(M());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << 3*m << endl;
    output << "property double x" << endl; 
    output << "property double y" << endl; 
    output << "property double z" << endl;
    output << "property double w" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<m;++i) {
        for(size_t j(0);j<3;++j){
            size_t k = triangles[i][j];
            Vertex<double> vec(vertices[k].x,vertices[k].y,vertices[k].z,values[3*i+j]);
            output.write((char*) (&vec),sizeof(Vertex<double>));
        }
    }
    for(size_t i(0);i<m;++i){
        uint_vec vec(3*i,3*i+1,3*i+2);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}



struct float_4 {
    float x,y,z,w;
};

struct double_4 {
    double x,y,z,w;
};

struct float_3 {
    float x,y,z;
};

struct uint_3 {
    unsigned int i,j,k;
};

struct int_3 {
    int i,j,k;
};

enum Property {
    FLOAT,
    DOUBLE,
    UINT,
    INT,
    CHAR,
    LIST_CHAR_UINT,
    LIST_CHAR_INT
};

#include <array>

array<size_t,7> propSize {
    sizeof(float),
    sizeof(double),
    sizeof(unsigned int),
    sizeof(int),
    sizeof(char),
    sizeof(unsigned int),
    sizeof(int)
};


vector<Property> read_properties(ifstream& input){
    vector<Property> result;

    string line;
    size_t pos(input.tellg());
    while(getline(input,line)) {
        //cout << "test" << endl;
        stringstream sstream(line);
        string word;
        sstream >> word;

        if(word == "element" or word == "end_header") {
            input.seekg(pos);
            return result;
        }

        else if(word == "property") {
            string type;
            //cout << "property   ";
            sstream >> type;

            if(type == "float" or type == "float32") {
                result.push_back(FLOAT);
                //cout << "float" << endl;
            }

            else if(type == "double" or type == "float64") {
                result.push_back(DOUBLE);
                //cout << "double" << endl;
            }

            else if(type == "uchar" or type == "uint8" or type == "char") {
                result.push_back(CHAR);
                //cout << "char" << endl;
            }

            else if(type == "list") {
                sstream >> type;
                if(type == "uchar" or type == "uint8") {
                    sstream >> type;
                    if(type == "uint" or type == "uint32") {
                        result.push_back(LIST_CHAR_UINT);
                        //cout << "list char uint" << endl;
                    }
                    else if(type == "int" or type == "int32") {
                        result.push_back(LIST_CHAR_INT);
                        //cout << "list char int" << endl;
                    }
                } else {
                    cerr << "unsupported type for polygon size: '" << type << "' instead of 'uchar' or 'uint8" << endl;
                }
                
            }
        }

        else {
            input.seekg(pos);
            cerr << "Invalid description string found (not property)" << endl;
            return result;
        }

        pos = input.tellg();
    }
    cerr << "couldn't continue reading." << endl;
    return result;
}


void Mesh::import_ply(string filename,vector<real>& values) {

    // first, clean everything up.
    vertices.clear();
    triangles.clear();
    values.clear();


    ifstream input(filename);
    string line;
    bool reading(true);

    size_t num_vertex(0);
    size_t num_faces(0);

    vector<Property> vertexprops;
    vector<Property> faceprops;


    while(reading and getline(input,line)) {
        if(line.find("end_header") != string::npos) { // meaning if end_header is found
            reading = false;
            //cout << "finished reading header." << endl;
        } else {
            string word;
            stringstream sstream(line);

            sstream >> word;
            if(word == "element") {
                sstream >> word;
                
                if(word == "vertex"){
                    //cout << "vertex   " << endl;
                    sstream >> num_vertex;
                    vertexprops = read_properties(input);
                }

                else if(word == "face"){
                    //cout << "face   " << endl;
                    sstream >> num_faces;
                    faceprops = read_properties(input);
                }

                else {
                    cerr << "not known";
                }
            }
        }
    }

    size_t pos(input.tellg());

    input.close();


    // reading the binary part: vertices

    ifstream binarystream(filename,ifstream::binary);
    if(binarystream){

        vector<char> text(pos);
        binarystream.read(text.data(),text.size());

#ifdef VERBOSE
        for(char elm : text) {
            cout << elm;
        }
        cout << "-------------------" << endl;
#endif

        binarystream.seekg(pos);

        size_t vertex_size(0);
        for(Property prop : vertexprops) {
            vertex_size += propSize[prop];
        }

        if(vertexprops[0] == FLOAT and vertexprops[1] == FLOAT and vertexprops[2] == FLOAT){
            if(vertexprops[3] == FLOAT){
#ifdef VERBOSE
                cout << "reading 4 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_4 result = parse_binary<float_4>(binarystream,vertex_size);
                    vertices.push_back(vec3(result.x,result.y,result.z));
                    values.push_back(result.w);
                }
            }
            else {
#ifdef VERBOSE
                cout << "reading 3 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_3 result = parse_binary<float_3>(binarystream,vertex_size);
                    vertices.push_back(vec3(result.x,result.y,result.z));
                }
            }
        }
        else if(vertexprops[0] == DOUBLE and vertexprops[1] == DOUBLE and vertexprops[2] == DOUBLE and vertexprops[3] == DOUBLE){
#ifdef VERBOSE
            cout << "reading 4 doubles" << endl;
#endif
            for(size_t i(0);i<num_vertex;++i){
                double_4 result = parse_binary<double_4>(binarystream,vertex_size);
                vertices.push_back(vec3(result.x,result.y,result.z));
                values.push_back(result.w);
            }
        }
        else {
            cerr << "number of detected floats/doubles not compatible with any known pattern." << endl;
            throw;
        }

        // Now reading the triangle indices:

        size_t face_rest_size(sizeof(uint_3));
        for(size_t i(1);i<faceprops.size();++i) {
            face_rest_size += propSize[faceprops[i]];
        }

        if(faceprops[0] == LIST_CHAR_UINT){
#ifdef VERBOSE
            cout << "reading 3 unsigned ints" << endl;
#endif
            for(size_t i(0);i<num_faces;++i){
                char num;
                binarystream.read(&num,1);
                if(num == 3) {
                    uint_3 result = parse_binary<uint_3>(binarystream,face_rest_size);
                    triangles.push_back(Triplet(result.i,result.j,result.k));
                }
                else {
                    cerr << "invalid number of indices: no polygons larger than triangles allowed!" << endl;
                    throw;
                }
            }
    
        }
        else if (faceprops[0] == LIST_CHAR_INT){
#ifdef VERBOSE
            cout << "reading 3 ints" << endl;
#endif
            for(size_t i(0);i<num_faces;++i){
                char num;
                binarystream.read(&num,1);
                if(num == 3) {
                    int_3 result = parse_binary<int_3>(binarystream,face_rest_size);
                    triangles.push_back(Triplet(result.i,result.j,result.k));
                }
                else {
                    cerr << "invalid number of indices: no polygons larger than triangles allowed!" << endl;
                    throw;
                }
            }
        }
        else {
            cerr << "Unknown pattern for face indices: first property of face element shall be an integer list beginning with a char" << endl;
            throw;
        }
        
        
        
    }
#ifdef VERBOSE
    cout << "finished reading '"+filename+"'." << endl;
#endif

}
