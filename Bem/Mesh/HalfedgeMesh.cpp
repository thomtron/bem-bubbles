#include "HalfedgeMesh.hpp"
#include "Tuplet.hpp"

#include <vector>
#include <algorithm> // for sort()
#include <numeric>   // for iota()

using namespace std;

namespace Bem {

const size_t HalfedgeMesh::npos = -1;

HalfedgeMesh::~HalfedgeMesh() {
    clear();
}

void HalfedgeMesh::clear() {
    for(Halfedge* elm : trigs) {
        delete elm->next->next;
        delete elm->next;
        delete elm;
    }
    for(Halfedge* elm : bounds) {
        Halfedge* u(elm->next);
        while(u != elm) {
            Halfedge* tmp = u->next;
            delete u;
            u = tmp;
        }
        delete elm;
    }
    verts.clear();
    edges.clear();
    trigs.clear();
    bounds.clear();
}

HalfedgeMesh::HalfedgeMesh(Mesh const& other){
    generate_halfedges(*this,other);
}

HalfedgeMesh::HalfedgeMesh(HalfedgeMesh const& other) {
    copy(other);
}
HalfedgeMesh& HalfedgeMesh::operator=(HalfedgeMesh const& other) {
    copy(other);
    return *this;
}

void HalfedgeMesh::copy(HalfedgeMesh const& other) {
    cout << "making copy." << endl;
    vpos = other.vpos;
    verts = vector<Halfedge*>(other.verts.size(),nullptr);
    edges = vector<Halfedge*>(other.edges.size(),nullptr);
    trigs = vector<Halfedge*>(other.trigs.size(),nullptr);
    for(Halfedge* elm : other.edges) {
        if(elm->trig == npos) { // if at boundary
            Halfedge* A = new Halfedge;

            A->vert = elm->vert;
            A->edge = elm->edge;
            A->trig = elm->trig;
            A->twin = A;

            verts[A->vert] = A;
            edges[A->edge] = A;
            trigs[A->trig] = A;
        } else {
            Halfedge* A = new Halfedge;
            Halfedge* B = new Halfedge;

            A->vert = elm->vert;
            A->edge = elm->edge;
            A->trig = elm->trig;
            A->twin = B;

            B->vert = elm->twin->vert;
            B->edge = elm->twin->edge;
            B->trig = elm->twin->trig;
            B->twin = A;

            verts[A->vert] = A;
            edges[A->edge] = A;
            trigs[A->trig] = A;

            verts[B->vert] = B;
            edges[B->edge] = B;
            trigs[B->trig] = B; 
        }
    }
    for(Halfedge* elm : other.trigs) {
        size_t trig_ind = elm->trig;
        Halfedge* A = edges[elm->edge];
        Halfedge* B = edges[elm->next->edge];
        Halfedge* C = edges[elm->next->next->edge];
        if(A->trig != trig_ind) A = A->twin;
        if(B->trig != trig_ind) B = B->twin;
        if(C->trig != trig_ind) C = C->twin;

        A->next = B;
        B->next = C;
        C->next = A;
    }
    for(Halfedge* elm : other.bounds) {
        Halfedge* half = edges[elm->edge];
        if(half->twin->trig == HalfedgeMesh::npos) half = half->twin;
        bounds.push_back(half);

        Halfedge* start = half->twin;
        Halfedge* curr = start;
        Halfedge* curr_before = half;
        do {
            curr = curr->next->twin;
            if(curr->trig == HalfedgeMesh::npos) {
                curr->next = curr_before;
                curr_before = curr;
                curr = curr->twin;
            } 
        } while(curr != start);
    }
}

HalfedgeMesh generate_halfedges(Mesh const& mesh) {
    HalfedgeMesh result;
    generate_halfedges(result,mesh);
    return result;
}

void generate_halfedges(HalfedgeMesh& result, Mesh const& mesh) {
    result.clear();
    
    // output arrays
    result.vpos = mesh.verts;
    vector<Halfedge*>& verts = result.verts;
    verts = vector<Halfedge*>(mesh.verts.size(),nullptr);
    vector<Halfedge*>& trigs = result.trigs;

    vector<Tuplet> edges;
    vector<Halfedge*> edge_halfs;

    for(size_t i(0);i<mesh.trigs.size();++i) {
        Triplet t(mesh.trigs[i]);

        Halfedge* A = new Halfedge;
        Halfedge* B = new Halfedge;
        Halfedge* C = new Halfedge;

        // here all relations between faces,vertices and halfedges are created

        A->next = B;
        A->vert = t.a;
        A->trig = i;

        B->next = C;
        B->vert = t.b;
        B->trig = i;

        C->next = A;
        C->vert = t.c;
        C->trig = i;

        trigs.push_back(A);

        verts[t.a] = A; // these are overwritten multiple times in this loop
        verts[t.b] = B;
        verts[t.c] = C;

        // the "glue" between the triangles is given by the edges, for which 
        // we have to generate a list using the Tuplet class

        edges.push_back(Tuplet(t.a,t.b));
        edges.push_back(Tuplet(t.b,t.c));
        edges.push_back(Tuplet(t.c,t.a));

        edge_halfs.push_back(A);
        edge_halfs.push_back(B);
        edge_halfs.push_back(C);
    }

    size_t K(edges.size());
    vector<size_t> edge_inds(K);
    iota(edge_inds.begin(),edge_inds.end(),0);
    sort(edge_inds.begin(),edge_inds.end(),[&edges](size_t a,size_t b) { return edges[a] < edges[b]; });

    size_t edge_ind(0);
    for(size_t i(0);i<K;++i) {

        if(i == K-1 or edges[edge_inds[i]]<edges[edge_inds[i+1]]) { 
            // if the edge was pushed back two times (from two triangles), the edge is
            // an interior edge. The present condition evaluates false, since the two 
            // indices are identical. In the other case, the edge was pushed back only
            // once and therefore we have a boundary edge, the case treated here:

            size_t ind  = edge_inds[i];
            Halfedge* A = edge_halfs[ind];

            A->twin = A;
            A->edge = edge_ind;

            result.edges.push_back(A);

        } else {

            size_t ind1 = edge_inds[i];
            i++;
            size_t ind2 = edge_inds[i];

            Halfedge* A = edge_halfs[ind1];
            Halfedge* B = edge_halfs[ind2];

            // now we can add the final informations to the halfedges

            A->twin = B;
            B->twin = A;

            A->edge = edge_ind;
            B->edge = edge_ind;

            result.edges.push_back(A);
        }

        edge_ind++;

    }

    for(Halfedge* u : result.verts) {
        Halfedge* start(u);
        Halfedge* first_boundary(nullptr);
        Halfedge* last_boundary(nullptr);

        //size_t k(0);
        do {
            if(u == u->twin) { // handle boundary edges! (add corresponding halfedges)
                //cout << "  > edge-piece nr. " << k++ << endl;
                Halfedge* B = new Halfedge;
                B->twin = u;
                B->trig = HalfedgeMesh::npos;
                B->next = last_boundary;
                B->edge = u->edge;
                B->vert = u->next->vert;

                u->twin = B;
                u = u->next;

                last_boundary = B;
                if(first_boundary == nullptr) {
                    first_boundary = B;
                    // since this part of the code is called once for each boundary,
                    // we can append now an element to the bounds array
                    result.bounds.push_back(B);
                } 

            } else {
                u = u->twin->next;
            }
        } while(u != start);
        if(first_boundary != nullptr) {
            first_boundary->next = last_boundary;
        }
    }
}

Mesh generate_mesh(HalfedgeMesh const& mesh) {
    Mesh result;
    generate_mesh(result,mesh);
    return result;
}

void generate_mesh(Mesh& result, HalfedgeMesh const& mesh) {
    result.clear();
    result.verts = mesh.vpos;
    for(Halfedge* elm : mesh.trigs) {
        Triplet t(elm->vert,elm->next->vert,elm->next->next->vert);
        result.trigs.push_back(t);
    }
}

bool at_boundary_u(Halfedge* vert) {
    Halfedge* u(vert);
    do {
        if(u->trig == HalfedgeMesh::npos or u->twin->trig == HalfedgeMesh::npos) return true;
        u = u->twin->next;
    } while(u != vert);
    return false;
}

// a function to check the validity of a HalfedgeMesh (check that all pointers 
// point to an element and that they obey the connectivity rules). It was mainly
// introduced for debugging.
bool HalfedgeMesh::check_validity() const {
#ifdef VERBOSE
    cout << "checking validity of generated mesh: ";
#endif
    // number of errors
    size_t error_nnn(0);
    size_t error_trigind(0);
    size_t error_tt(0);
    size_t error_edgeind(0);
    size_t error_vertind(0);
    size_t error_valence(0);

    //cout << "a" << endl;

    for(Halfedge* elm : trigs) {
        if(elm->next->next->next != elm) error_nnn++;           // check that next-pointers connect correctly
        if(elm->next->trig != elm->trig) error_trigind++;       // check that halfedges point to same triangle
        if(elm->next->next->trig != elm->trig) error_trigind++; // -^
    }
    //cout << "b" << endl;
    for(Halfedge* elm : edges) {
        if(elm->twin->twin != elm) error_tt++;                  // check that twin of twin is element itself
        if(edges[elm->edge] != elm) error_edgeind++;       // check that edge indices are correct
        if(edges[elm->twin->edge] != elm) error_edgeind++; // check that edge indices are correct
    }
    //cout << "c" << endl;


    // the following 18 lines are for checking that the valence number computed by two different methods are consistent
    vector<size_t> valences(verts.size(),0);

    for(Halfedge* elm : trigs) {
        valences[elm->vert]++;
        valences[elm->next->vert]++;
        valences[elm->next->next->vert]++;
    }
    //cout << "d" << endl;

    // add one for all vertices at boundaries
    for(Halfedge* elm : bounds) {
        Halfedge* u(elm);
        do {
            valences[u->vert]++;
            u = u->next; // slide along boundary
        } while(u != elm);
    }
    //cout << "e" << endl;

    size_t i(0);
    for(Halfedge* elm : verts) {
        i++;
        Halfedge* u(elm);
        size_t valence = 0;
        do {
            valence++;
            if(u->vert != elm->vert) {
                error_vertind++;
                if(at_boundary_u(verts[u->vert])) cout << "boundary: " << i << "/" << verts.size() << endl;
            } 
            u = u->twin->next;
        } while(u != elm);
        if(valence != valences[elm->vert]){
            error_valence++; 
            cout << valence << " -%- " << valences[elm->vert] << endl;
        } 
    }
    //cout << "f" << endl;


//#ifdef VERBOSE
    // print found errors and their occurence
    if(error_nnn>0)     cout << "error: next->next->next not unity  (" << error_nnn     << ")" << endl;
    if(error_trigind>0) cout << "error: bad triangle index          (" << error_trigind << ")" << endl;
    if(error_tt>0)      cout << "error: twin->twin not unity        (" << error_tt      << ")" << endl;
    if(error_edgeind>0) cout << "error: bad edge index              (" << error_edgeind << ")" << endl;
    if(error_valence>0) cout << "error: bad valence number          (" << error_valence << ")" << endl;
    if(error_vertind>0) cout << "error: bad vert index              (" << error_vertind << ")" << endl;
//#endif

    return (error_nnn     == 0 and
            error_trigind == 0 and
            error_edgeind == 0 and
            error_valence == 0 and
            error_vertind == 0 );
}


vector<vec3> generate_vertex_normals(HalfedgeMesh const& mesh) {
    vector<vec3> vertex_normals;

    size_t n(mesh.verts.size());

    for(size_t i(0);i<n;++i) {
        vec3 normal;
        Halfedge* u(mesh.verts[i]);
        Halfedge* v(u);
        do {
            vec3 B(mesh.vpos[u->next->vert]-mesh.vpos[u->vert]);
            vec3 C(mesh.vpos[u->next->next->vert]-mesh.vpos[u->vert]);
            normal += B.vec(C)*(1.0/(B.norm2()*C.norm2()));


            u = u->twin->next;
        } while(u != v);
        normal.normalize();
        vertex_normals.push_back(normal);
    }
    return vertex_normals;
}


} // namespace Bem