#include "HalfedgeMesh.hpp"
#include "Tuplet.hpp"

#include <vector>
#include <list>
#include <algorithm> // for sort()
#include <numeric>   // for iota()

using namespace std;

namespace Bem {

HalfedgeMesh::~HalfedgeMesh() {
    release();
}

void HalfedgeMesh::release() {
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

void HalfedgeMesh::copy(HalfedgeMesh const& other) { // can be optimized very probably

    generate_halfedges(*this,generate_mesh(other));
    /*
    cout << "making copy." << endl;

    for(Halfedge* elm : other.edges) {
        Halfedge* A = new Halfedge;
        Halfedge* B = new Halfedge;

        edges.push_back(A);

        A->edge = &edges.back();
        A->twin = B;

        B->edge = &edges.back();
        B->twin = A; 
    }
    for(Vertex const& elm : other.verts) {
        verts.push_back(elm);
        Vertex* vert(&verts.back());
        vert->half = 
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
    */
}

/*
HalfedgeMesh generate_halfedges(Mesh const& mesh) {

    HalfedgeMesh result;
    
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

            size_t ind = edge_inds[i];
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

    return result;
}
*/

HalfedgeMesh generate_halfedges(Mesh const& mesh) {
    HalfedgeMesh result;
    generate_halfedges(result,mesh);
    return result;
}

void generate_halfedges(HalfedgeMesh& result, Mesh const& mesh) {
    result.release();
    // output arrays

    vector<Vertex> verts;
    
    for(vec3 const& vec : mesh.verts) {
        verts.push_back({nullptr,vec,0});
    }

    //result.trigs = vector<Halfedge*>(mesh.trigs.size());

    
    vector<Tuplet> edges;
    vector<Halfedge*> edge_halfs;

    for(size_t i(0);i<mesh.trigs.size();++i) {
        Triplet t(mesh.trigs[i]);

        Halfedge* A = new Halfedge;
        Halfedge* B = new Halfedge;
        Halfedge* C = new Halfedge;

        // here all relations between faces,vertices and halfedges are created

        //result.trigs[i] = A;
        result.trigs.push_back(A);

        //result.verts[t.a].half = A; // these are overwritten multiple times in this loop
        //result.verts[t.b].half = B;
        //result.verts[t.c].half = C;

        verts[t.a].half = A; // these are overwritten multiple times in this loop
        verts[t.b].half = B;
        verts[t.c].half = C;

        A->next = B;
        A->trig = &result.trigs.back();
        A->edge = nullptr;
        A->twin = nullptr;

        B->next = C;
        B->trig = &result.trigs.back();
        B->edge = nullptr;
        B->twin = nullptr;

        C->next = A;
        C->trig = &result.trigs.back();
        C->edge = nullptr;
        C->twin = nullptr;

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

    for(size_t i(0);i<K;++i) {

        if(i == K-1 or edges[edge_inds[i]]<edges[edge_inds[i+1]]) { 
            // if the edge was pushed back two times (from two triangles), the edge is
            // an interior edge. The present condition evaluates false, since the two 
            // indices are identical. In the other case, the edge was pushed back only
            // once and therefore we have a boundary edge, the case treated here:

            size_t ind  = edge_inds[i];
            Halfedge* A = edge_halfs[ind];

            result.edges.push_back(A);

            A->twin = A;
            A->edge = &result.edges.back();

        } else {

            size_t ind1 = edge_inds[i];
            i++;
            size_t ind2 = edge_inds[i];

            Halfedge* A = edge_halfs[ind1];
            Halfedge* B = edge_halfs[ind2];

            // now we can add the final informations to the halfedges

            A->twin = B;
            B->twin = A;

            result.edges.push_back(A);

            
            A->edge = &result.edges.back();
            B->edge = &result.edges.back();

            //Halfedge** h = result.edges[0]->edge;

            //cout << A->edge << " - " << h << endl;

            //cout << result.get_index(result.edges,h) << ' ';// << endl;


        }
    }


    result.verts = list<Vertex>(verts.begin(),verts.end());
    
    for(Vertex& elm : result.verts) {
        Halfedge* u(elm.half);
        Halfedge* first_boundary(nullptr);
        Halfedge* last_boundary(nullptr);

        u = u->twin->next;
        //size_t k(0);
        while(u != elm.half) {
            if(u == u->twin) { // handle boundary edges! (add corresponding halfedges)
                //cout << "  > edge-piece nr. " << k++ << endl;
                Halfedge* B = new Halfedge;
                B->twin = u;
                B->trig = nullptr;
                B->next = last_boundary;
                B->edge = u->edge;

                u->twin = B;
                u = u->next;

                last_boundary = B;
                if(first_boundary == nullptr) first_boundary = B;

            } else {
                u = u->twin->next;
            }
        }
        if(first_boundary != nullptr) {
            first_boundary->next = last_boundary;
        }
    }

    for(Vertex& elm : result.verts) {
        Halfedge* u(elm.half);
        u->vert = &elm;
        u = u->twin->next;
        while(u != elm.half) {
            u->vert = &elm;
            u = u->twin->next;
        }
    }

    result.update_vertex_indices();
}

Mesh generate_mesh(HalfedgeMesh const& mesh) {
    Mesh result;
    generate_mesh(result,mesh);
    return result;
}

void generate_mesh(Mesh& result, HalfedgeMesh const& mesh) {

    for(Vertex const& elm : mesh.verts)
        result.verts.push_back(elm.pos);

    for(Halfedge* elm : mesh.trigs) {
        Triplet t;

        //vector<Vertex>::iterator const& it(*(elm->vert));
        //auto index = (elm->vert - &mesh.verts[0])/sizeof(Vertex); // mesh.verts.begin() + 
        //auto index = distance(mesh.verts.begin(),vector<Vertex>::iterator(&mesh.verts[7])); //vector<const Vertex>::iterator(elm->vert));

        t.a = elm->vert->index;
        t.b = elm->next->vert->index;
        t.c = elm->next->next->vert->index;

        result.trigs.push_back(t);
    }
}

void HalfedgeMesh::update_vertex_indices() {
    size_t i(0);
    for(Vertex& elm : verts) {
        elm.index = i;
        i++;
    }
}


} // namespace Bem