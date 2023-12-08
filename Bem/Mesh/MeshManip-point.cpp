#include "MeshManip.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <omp.h> // for project

#include "MeshIO.hpp"
#include "FittingTool.hpp"

using namespace std;

namespace Bem {

// a function to check the validity of a HalfedgeMesh (check that all pointers 
// point to an element and that they obey the connectivity rules). It was mainly
// introduced for debugging.
void check_validity(HalfedgeMesh const& mesh) {
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

    for(Halfedge* elm : mesh.trigs) {
        if(elm->next->next->next != elm) error_nnn++;           // check that next-pointers connect correctly
        if(elm->next->trig != elm->trig) error_trigind++;       // check that halfedges point to same triangle
        if(elm->next->next->trig != elm->trig) error_trigind++; // -^
    }
    for(Halfedge* const& elm : mesh.edges) {
        //cout << mesh.get_index(mesh.edges,&elm) << endl;
        if(elm->twin->twin != elm) error_tt++;                   // check that twin of twin is element itself
        if(elm->edge != elm->twin->edge) error_edgeind++;        // check that edge indices are correct
        if(*(elm->edge) != elm and *(elm->edge) != elm->twin) {
            error_edgeind++; // check that edge indices are correct
            if(elm->trig == nullptr) cout << "elm->trig is nullptr" << endl;
            if(elm->twin->trig == nullptr) cout << "elm->twin->trig is nullptr" << endl;
            if(elm->edge == nullptr) cout << "nullptr" << endl;
            //cout << mesh.get_index(mesh.edges,&elm) << " -%- " << mesh.get_index(mesh.edges,elm->edge)<< endl;
            
        }
    }

    /*
    // the following 18 lines are for checking that the valence number computed by two different methods are consistent
    vector<size_t> valences(mesh.verts.size(),0);
    
    for(Halfedge* elm : mesh.trigs) {
        valences[mesh.get_index(mesh.verts,elm->vert)]++;
        valences[mesh.get_index(mesh.verts,elm->next->vert)]++;
        valences[mesh.get_index(mesh.verts,elm->next->next->vert)]++;
    }

    for(Vertex const& elm : mesh.verts) {
        Halfedge* half(elm.half);
        Halfedge* u(half);
        size_t valence = 0;
        do {
            valence++;
            if(u->vert != half->vert) error_vertind++;
            u = u->twin->next;
        } while(u != half);
        if(valence != valences[mesh.get_index(mesh.verts,half->vert)]) {
            error_valence++; //cout << "v-e: " << valence << ", exp: " << valences[elm->vert] << endl; //error_valence++;
            //cout << int(valence) - int(valences[mesh.get_index(mesh.verts,half->vert)]) << endl;
        } 
    }*/

    // print found errors and their occurence

#ifdef VERBOSE
    if(error_nnn     == 0 and
       error_trigind == 0 and
       error_edgeind == 0 and
       error_valence == 0 and
       error_vertind == 0 ) cout << "okay";

    cout << endl;
#endif

    if(error_nnn>0)     cout << "error: next->next->next not unity  (" << error_nnn     << ")" << endl;
    if(error_trigind>0) cout << "error: bad triangle index          (" << error_trigind << ")" << endl;
    if(error_tt>0)      cout << "error: twin->twin not unity        (" << error_tt      << ")" << endl;
    if(error_edgeind>0) cout << "error: bad edge index              (" << error_edgeind << ")" << endl;
    if(error_valence>0) cout << "error: bad valence number          (" << error_valence << ")" << endl;
    if(error_vertind>0) cout << "error: bad vert index              (" << error_vertind << ")" << endl;
}


void split_edge(HalfedgeMesh& mesh,Halfedge* edge) {
    Halfedge* edge01(edge);
    Halfedge* edge10(edge01->twin);

    if(edge01->trig == nullptr or edge10->trig == nullptr) { // if at edge
        // a remplir
        //cout << "edge" << endl;
    } else {

        //mesh.vpos.push_back(0.5*(edge01->vert->pos+edge01->next->vert->pos));
        //mesh.verts.push_back(edge10);

        mesh.verts.push_back({edge10,0.5*(edge01->vert->pos+edge01->next->vert->pos),mesh.verts.size()});
        //size_t new_index(mesh.verts.size()-1);

        //Halfedge** new_vert() // damit, push_back invalidates again pointers!

        // S for straight, A, B for the two sides 
        // _0 are Halfedges going out from new_index
        //size_t K(mesh.edges.size());
        //size_t M(mesh.trigs.size());

        Vertex* new_vert = &mesh.verts.back();

        Halfedge* S_0 = new Halfedge;
        Halfedge* S_1 = new Halfedge;
        Halfedge* A_0 = new Halfedge;
        Halfedge* A_1 = new Halfedge;
        Halfedge* B_0 = new Halfedge;
        Halfedge* B_1 = new Halfedge;

        mesh.edges.push_back(S_0);
        Halfedge** edge_S = &mesh.edges.back();
        mesh.edges.push_back(A_0);
        Halfedge** edge_A = &mesh.edges.back();
        mesh.edges.push_back(B_0);
        Halfedge** edge_B = &mesh.edges.back();

        mesh.trigs.push_back(S_0);
        Halfedge** trig_A = &mesh.trigs.back();
        mesh.trigs.push_back(S_1);
        Halfedge** trig_B = &mesh.trigs.back();

        S_0->edge = edge_S;
        S_0->trig = trig_A;
        S_0->next = edge01->next;
        S_0->vert = new_vert;
        S_0->twin = S_1;

        S_1->edge = edge_S;
        S_1->trig = trig_B;
        S_1->next = B_0;
        S_1->vert = edge01->next->vert;
        S_1->twin = S_0;
    
        A_0->edge = edge_A;
        A_0->trig = edge01->trig;
        A_0->next = edge01->next->next;
        A_0->vert = new_vert;
        A_0->twin = A_1;
        
        A_1->edge = edge_A;
        A_1->trig = trig_A;
        A_1->next = S_0;
        A_1->vert = edge01->next->next->vert;
        A_1->twin = A_0;

        B_0->edge = edge_B;
        B_0->trig = trig_B;
        B_0->next = edge10->next->next;
        B_0->vert = new_vert;
        B_0->twin = B_1;

        B_1->edge = edge_B;
        B_1->trig = edge10->trig;
        B_1->next = edge10;
        B_1->vert = edge10->next->next->vert;
        B_1->twin = B_0;


        edge01->next->trig = trig_A;
        edge01->next->next = A_1;
        edge01->next = A_0;

        edge10->next->next->trig = trig_B;
        edge10->next->next->next = S_1;
        edge10->next->next = B_1;
        edge10->vert = new_vert;

        S_1->vert->half = S_1;
        *(edge01->trig) = edge01;
        *(edge10->trig) = edge10;
    }
}


// split_edges splits an edge if it is longer than L_max into two smaller edges, creating a new vertex
// at the midpoint of the original edge. Another set of two edges and two triangles have to be added too
// to render the triange mesh valid (without holes)
void split_edges    (HalfedgeMesh& mesh, real L_max) {
    check_validity(mesh);
#ifdef VERBOSE
    cout << "SPLIT-EDGES" << endl;
    size_t num_t_init = mesh.trigs.size();
#endif
    
    // we work with square lengths since we only compare lengths qualitatively and 
    // the computation of the square root can be prevented like that.
    vector<real> lengths;
    real L_max_2 = L_max*L_max;
    
    for(Halfedge* halfedge : mesh.edges) {
        lengths.push_back((halfedge->next->vert->pos-halfedge->vert->pos).norm2());
    }

    // reorder the halfedge mesh and lengths such that the longest edges come first
    // create an index list
    vector<size_t> indices(lengths.size());
    iota(indices.begin(),indices.end(),0);
    // sort the index list according to the lengths of the edges
    sort(indices.begin(),indices.end(),[&lengths](size_t a,size_t b) { return lengths[a] > lengths[b]; });

    vector<Halfedge*> edges_tmp(mesh.edges.begin(),mesh.edges.end());
    vector<real> lengths_tmp = lengths;

    // apply the permutation stored in indices to lengths and the halfedgemesh
    size_t i(0);
    for(Halfedge*& half : mesh.edges) {
        size_t k = indices[i];
        half = edges_tmp[k];
        half->edge = &half;
        half->twin->edge = &half;
        lengths[i] = lengths_tmp[k];
        i++;
    }  

    // we want to check each edge only once
    // (the newly generated edges will be appended at the end and we won't check them in this pass)
    list<Halfedge*>::iterator edge = mesh.edges.begin();
    size_t I(lengths.size());
    for(size_t i(0);i<I;++i) {
        
        // check whether edge is longer than L_max
        if(lengths[i] > L_max_2)
            split_edge(mesh,*edge);

        ++edge;
    }
#ifdef VERBOSE
    cout << endl << "enhanced mesh from " << num_t_init << " to " << mesh.trigs.size() << " triangles." << endl;
#endif
}


// same function as above but this time with an array 'curvature' that 
// controls the maximum edgelength at each vertex the multiplicator gives a global
// factor on this maximum length. the curvature values are inversely proportional
// to the maximum length. Alternative: just pass a vector<real>& max_length as argument
// and leave the conversion of curvature to length outside. the multiplicator is 
// useful altough!
void split_edges    (HalfedgeMesh& mesh, vector<real>& max_edgelength) {
    assert(max_edgelength.size() == mesh.verts.size());
    check_validity(mesh);
#ifdef VERBOSE
    cout << "SPLIT-EDGES" << endl;
    size_t num_t_init = mesh.trigs.size();
#endif
    // we work with square lengths since we only compare lengths qualitatively and 
    // the computation of the square root can be prevented like that.
    vector<real> lengths;
    
    for(Halfedge* halfedge : mesh.edges) {
        lengths.push_back((halfedge->next->vert->pos-halfedge->vert->pos).norm2());
    }

    // reorder the halfedge mesh and lengths such that the longest edges come first
    // create an index list
    vector<size_t> indices(lengths.size());
    iota(indices.begin(),indices.end(),0);
    // sort the index list according to the lengths of the edges
    sort(indices.begin(),indices.end(),[&lengths](size_t a,size_t b) { return lengths[a] > lengths[b]; });

    vector<Halfedge*> edges_tmp(mesh.edges.begin(),mesh.edges.end());
    vector<real> lengths_tmp = lengths;
    vector<real> maxlengths_tmp = max_edgelength;

    // apply the permutation stored in indices to lengths and the halfedgemesh
    size_t i(0);
    for(Halfedge*& half : mesh.edges) {
        size_t k = indices[i];
        half = edges_tmp[k];
        half->edge = &half;
        half->twin->edge = &half;
        lengths[i] = lengths_tmp[k];
        i++;
    }
    
    mesh.update_vertex_indices();
    // we want to check each edge only once
    // (the newly generated edges will be appended at the end and we won't check them in this pass)
    list<Halfedge*>::iterator edge = mesh.edges.begin();
    size_t I(lengths.size());
    for(size_t i(0);i<I;++i) {
        // interpolate the max_edgelength
        real L_max = 0.5*(max_edgelength[(*edge)->vert->index] + max_edgelength[(*edge)->twin->vert->index]);
        real L_max_2 = L_max*L_max; // make square
        // check whether edge is longer than L_max
        if(lengths[i] > L_max_2) {
            split_edge(mesh,*edge);
            max_edgelength.push_back(L_max);
        }

        ++edge; // we can do this, since edge array is only growing
    }

#ifdef VERBOSE
    cout << endl << "enhanced mesh from " << num_t_init << " to " << mesh.trigs.size() << " triangles." << endl;
#endif
}

// for the collapse_edge function we need to be able to remove individual elements
// from the HalfedgeMesh. Here are some functions provided that render this task a bit easier

void remove_edge(vector<Halfedge*>& edges,size_t edge_ind) {
    /*
    if(edge_ind != edges.size()-1) {
        swap(edges[edge_ind],edges.back());
        Halfedge* u(edges[edge_ind]);
        u->edge = edge_ind;
        u->twin->edge = edge_ind;
    }
    edges.pop_back();
    */
}

void remove_vertex(vector<Halfedge*>& verts,size_t vert_ind) {
    /*
    if(vert_ind != verts.size()-1) {
        swap(verts[vert_ind],verts.back());
        Halfedge* u(verts[vert_ind]);
        Halfedge* v(u);
        do {
            v->vert = vert_ind;
            v = v->twin->next;
        } while(v != u);
    }
    verts.pop_back();
    */
}

void remove_face(vector<Halfedge*>& faces,size_t face_ind) {
    /*
    if(face_ind != faces.size()-1) {
        swap(faces[face_ind],faces.back());
        Halfedge* u(faces[face_ind]);
        u->trig = face_ind;
        u->next->trig = face_ind;
        u->next->next->trig = face_ind;
    }
    faces.pop_back();
    */
}

// the collapse_edges function tests whether an edge is shorter than L_min and, 
// in this case, remove this edge by merging its two corresponding vertices into one.
// Therefore its adjacent triangles and two of the edges that would become doubled
// have to be removed too.
void collapse_edges (HalfedgeMesh& mesh, real L_min) {
    check_validity(mesh);
#ifdef VERBOSE
    cout << "COLLAPSE-EDGES" << endl;
    size_t num_t_init = mesh.trigs.size();
#endif

    // we work with square lengths since we only compare lengths qualitatively and 
    // the computation of the square root can be prevented like that.
    vector<real> lengths;
    real L_min_2 = L_min*L_min;
    
    for(Halfedge* halfedge : mesh.edges) {
        lengths.push_back((halfedge->next->vert->pos-halfedge->vert->pos).norm2());
    }

    // reorder the halfedge mesh and lengths such that the longest edges come first
    // create an index list
    vector<size_t> indices(lengths.size());
    iota(indices.begin(),indices.end(),0);
    // sort the index list according to the lengths of the edges
    sort(indices.begin(),indices.end(),[&lengths](size_t a,size_t b) { return lengths[a] < lengths[b]; });

    vector<Halfedge*> edges_tmp(mesh.edges.begin(),mesh.edges.end());
    vector<real> lengths_tmp = lengths;

    // apply the permutation stored in indices to lengths and the halfedgemesh
    size_t i(0);
    for(Halfedge*& half : mesh.edges) {
        size_t k = indices[i];
        half = edges_tmp[k];
        half->edge = &half;
        half->twin->edge = &half;
        lengths[i] = lengths_tmp[k];
        i++;
    }  

    // looping through the edges while the edges-array becomes smaller
    //size_t k(0);
    for(list<Halfedge*>::iterator edgeit(mesh.edges.begin());edgeit != mesh.edges.end();++edgeit) {

        Halfedge* edge(*edgeit);
        //cout << edge << endl;
        //cout << "okay..." << endl;

        // check whether the current edge is shorter than the minimum length
        if((edge->next->vert->pos-edge->vert->pos).norm2() < L_min_2) {

            Halfedge* halfedge_A = edge;
            Halfedge* halfedge_B = halfedge_A->twin;

            if(halfedge_A->trig == nullptr or halfedge_B->trig == nullptr) {

            } else {

                //cout << k++ << endl;
                Vertex* vert_0 = halfedge_A->vert;
                Vertex* vert_1 = halfedge_B->vert;

                Vertex* vert_A = halfedge_A->next->next->vert;
                Vertex* vert_B = halfedge_B->next->next->vert;

                size_t num_paths(0);

                Halfedge* finder = vert_0->half;
                Halfedge* init = finder;
                do {
                    Halfedge* second = finder->next;
                    Halfedge* second_init = second;
                    do {
                        if(second->next->vert == vert_1) {
                            num_paths++;
                        }
                        second = second->twin->next;
                    } while(second != second_init);

                    finder = finder->twin->next;
                } while(finder != init);


                // Ensure that we wont produce "triangle appendices".
                // The edge is only collapsed, if it is geometrically admissible. If 
                // a triangle gets turned over or an appendix of zero
                // thickness would be produced (which is prevented by ensuring that
                // there are exactly two paths composed of two edges connecting the 
                // two vertices vert_0 and vert_1), we won't collapse the edge.
                // For this reason it is important to apply this function iteratively 
                // for good results.
                if(num_paths == 2) {

                    // test whether all triangles are still
                    // facing (approximately) in the right direction

                    vec3 new_pos(0.5*(vert_0->pos + vert_1->pos));
                    vec3 old_pos(vert_1->pos);

                    bool valid = true;

                    double limit_value = 0.8;

                    Halfedge* u(halfedge_A->next->twin);
                    do {
                        if(u == halfedge_B->next->next) {
                            u = halfedge_B->next->twin;
                            old_pos = vert_0->pos;
                        }
                        vec3 a = u->vert->pos - u->next->next->vert->pos;
                        vec3 b = u->next->next->vert->pos - new_pos;
                        vec3 c = u->next->next->vert->pos - old_pos;

                        c = c.vec(a);
                        c.normalize();
                        b = b.vec(a);
                        b.normalize();

                        if(b.dot(c)<limit_value) valid = false;
                        u = u->next->twin;

                    } while(valid and u != halfedge_A->next->next);


                    if(valid) {
                        //cout << "okay start" << endl;

                        // remove the elements

                        // remove the three edges
                        mesh.edges.erase(pointerToIter(halfedge_A->next->edge));
                        mesh.edges.erase(pointerToIter(halfedge_B->next->next->edge));
                        mesh.edges.erase(pointerToIter(halfedge_A->edge));


                        // remove the two faces
                        mesh.trigs.erase(pointerToIter(halfedge_A->trig));
                        mesh.trigs.erase(pointerToIter(halfedge_B->trig));


                        // make sure that vert_0 has an halfedge index that wont be removed
                        // the same for vert_A and vert_B
                        vert_0->half = halfedge_A->next->next->twin;
                        vert_A->half = halfedge_A->next->twin;
                        vert_B->half = halfedge_B->next->twin;

                        size_t index(vert_1->index); // !! brauchis? (für lengths-array [andere funktion])
                        mesh.verts.erase(pointerToIter(vert_1));

                        // adjust the real vertices of the mesh
                        vert_0->pos = new_pos;

                        // move pointers from vert_1 to vert_0
                        Halfedge* current = halfedge_A->next;
                        while( current != halfedge_B) {
                            current->vert = vert_0;
                            current = current->twin->next;
                        }

                        //cout << "okay 1" << endl;
                        // make sure that the edges of triangle A and B that WEREN'T removed, 
                        // will still have a valid halfedge pointer
                        *(halfedge_A->next->next->edge) = halfedge_A->next->next->twin;
                        *(halfedge_B->next->edge) = halfedge_B->next->twin;

                        //cout << "okay 2" << endl;
                        // now handle twins (before we destroy structure)
                        halfedge_A->next->next->twin->twin = halfedge_A->next->twin;
                        halfedge_A->next->twin->twin = halfedge_A->next->next->twin;
                        halfedge_A->next->twin->edge = halfedge_A->next->next->edge;

                        halfedge_B->next->next->twin->twin = halfedge_B->next->twin;
                        halfedge_B->next->twin->twin = halfedge_B->next->next->twin;
                        halfedge_B->next->next->twin->edge = halfedge_B->next->edge;

                        //cout << "okay 3" << endl;
                        // finally remove the halfedges

                        delete halfedge_A->next->next;
                        delete halfedge_A->next;
                        delete halfedge_A;

                        delete halfedge_B->next->next;
                        delete halfedge_B->next;
                        delete halfedge_B;

                        //cout << "okay 4" << endl;
                        
                        edgeit = mesh.edges.begin();

                        
                        //cout << "okay end" << endl;
                    }
                }
            }
        }
    }
#ifdef VERBOSE
    cout << endl << "simplified mesh from " << num_t_init << " to " << mesh.trigs.size() << " triangles." << endl;
#endif
}

// same function as above, but, as in the case of split_edges, with additional parameters 
// for adaptive refinement of the surface.
void collapse_edges (HalfedgeMesh& mesh, vector<real>& curvature, real multiplicator) {
    assert(curvature.size() == mesh.verts.size());
    check_validity(mesh);
#ifdef VERBOSE
    cout << "COLLAPSE-EDGES" << endl;
    size_t num_t_init = mesh.trigs.size();
#endif
    vector<real> lengths;
    /*
    for(auto halfedge :mesh.edges) {
        lengths.push_back((mesh.vpos[halfedge->next->vert]-mesh.vpos[halfedge->vert]).norm2());
    }

    // reorder the halfedge mesh and lengths such that shortest edges come first
    
    vector<size_t> indices(lengths.size());
    iota(indices.begin(),indices.end(),0);
    sort(indices.begin(),indices.end(),[&lengths](size_t a,size_t b) { return lengths[a] < lengths[b]; });

    vector<Halfedge*> edges_tmp = mesh.edges;
    vector<real> lengths_tmp = lengths;

    for(size_t i(0);i<lengths.size();++i) {
        size_t k = indices[i];
        mesh.edges[i] = edges_tmp[k];
        mesh.edges[i]->edge = i;
        mesh.edges[i]->twin->edge = i;
        lengths[i] = lengths_tmp[k];
        
    }


    for(size_t k(0);k<mesh.edges.size();++k) {

        Halfedge* halfedge_A = mesh.edges[k];
        Halfedge* halfedge_B = halfedge_A->twin;

        if(halfedge_A == halfedge_B) {

        } else {

            real curv_0 = curvature[halfedge_A->vert];
            real curv_1 = curvature[halfedge_B->vert];

            real L_min_2 = 0.0;
            bool zero_curv = false;
            if(curv_0+curv_1 != 0.0) {
                L_min_2 = multiplicator*2.0/(curv_0+curv_1);
                L_min_2 = L_min_2*L_min_2;
            }

            if(lengths[k] < L_min_2 or zero_curv) {

                Halfedge* halfedge_A = mesh.edges[k];
                Halfedge* halfedge_B = halfedge_A->twin;

                size_t vert_0 = halfedge_A->vert;
                size_t vert_1 = halfedge_B->vert;

                size_t vert_A = halfedge_A->next->next->vert;
                size_t vert_B = halfedge_B->next->next->vert;

                size_t num_paths(0);

                Halfedge* finder = mesh.verts[vert_0];
                Halfedge* init = finder;
                do {
                    Halfedge* second = finder->next;
                    Halfedge* second_init = second;
                    do {
                        if(second->next->vert == vert_1) {
                            num_paths++;
                        }
                        second = second->twin->next;
                    } while(second != second_init);

                    finder = finder->twin->next;
                } while(finder != init);


                if(num_paths == 2) {

                    vec3 new_pos(0.5*(mesh.vpos[vert_0] + mesh.vpos[vert_1]));
                    vec3 old_pos(mesh.vpos[vert_1]);

                    bool valid = true;

                    double limit_value = 0.8;

                    Halfedge* u(halfedge_A->next->twin);
                    do {
                        if(u == halfedge_B->next->next) {
                            u = halfedge_B->next->twin;
                            old_pos = mesh.vpos[vert_0];
                        }
                        vec3 a = mesh.vpos[u->vert] - mesh.vpos[u->next->next->vert];
                        vec3 b = mesh.vpos[u->next->next->vert] - new_pos;
                        vec3 c = mesh.vpos[u->next->next->vert] - old_pos;

                        c = c.vec(a);
                        c.normalize();
                        b = b.vec(a);
                        b.normalize();

                        if(b.dot(c)<limit_value) valid = false;
                        u = u->next->twin;

                    } while(valid and u != halfedge_A->next->next);


                    if(valid) {

                        // remove the elements

                        // remove the three edges
                        remove_edge(mesh.edges,halfedge_A->next->edge);
                        remove_edge(mesh.edges,halfedge_B->next->next->edge);
                        remove_edge(mesh.edges,halfedge_A->edge);

                        swap(lengths[halfedge_A->next->edge],lengths.back());
                        lengths.pop_back();
                        swap(lengths[halfedge_B->next->next->edge],lengths.back());
                        lengths.pop_back();
                        swap(lengths[halfedge_A->edge],lengths.back());
                        lengths.pop_back();


                        // remove the two faces
                        remove_face(mesh.trigs,halfedge_A->trig);
                        remove_face(mesh.trigs,halfedge_B->trig);


                        // make sure that vert_0 has an halfedge index that wont be removed
                        // the same for vert_A and vert_B
                        mesh.verts[vert_0] = halfedge_A->next->next->twin;
                        mesh.verts[vert_A] = halfedge_A->next->twin;
                        mesh.verts[vert_B] = halfedge_B->next->twin;

                        remove_vertex(mesh.verts,vert_1);

                        // adjust the real vertices of the mesh
                        mesh.vpos[vert_0] = new_pos;
                        swap(mesh.vpos[vert_1],mesh.vpos.back());
                        mesh.vpos.pop_back();

                        // and the curvature values
                        curvature[vert_0] = 0.5*(curv_0+curv_1);
                        swap(curvature[vert_1],curvature.back());
                        curvature.pop_back();

                        // in case that we just swapped vert_0, adjust it here (since we use it later on)
                        if(vert_0 == mesh.verts.size()) vert_0 = vert_1;

                        // move pointers from vert_1 to vert_0
                        Halfedge* current = halfedge_A->next;
                        while( current != halfedge_B) {
                            current->vert = vert_0;
                            current = current->twin->next;
                        }

                        // make sure that the edges of triangle A and B that WEREN'T removed, 
                        // will still have a valid halfedge pointer
                        mesh.edges[halfedge_A->next->next->edge] = halfedge_A->next->next->twin;
                        mesh.edges[halfedge_B->next->edge] = halfedge_B->next->twin;

                        // now handle twins (before we destroy structure)
                        halfedge_A->next->next->twin->twin = halfedge_A->next->twin;
                        halfedge_A->next->twin->twin = halfedge_A->next->next->twin;
                        halfedge_A->next->twin->edge = halfedge_A->next->next->edge;

                        halfedge_B->next->next->twin->twin = halfedge_B->next->twin;
                        halfedge_B->next->twin->twin = halfedge_B->next->next->twin;
                        halfedge_B->next->next->twin->edge = halfedge_B->next->edge;

                        // finally remove the halfedges

                        delete halfedge_A->next->next;
                        delete halfedge_A->next;
                        delete halfedge_A;

                        delete halfedge_B->next->next;
                        delete halfedge_B->next;
                        delete halfedge_B;

                        // adjust lengths

                        Halfedge* half_0 = mesh.verts[vert_0];
                        current = half_0;
                        do {
                            lengths[current->edge] = (mesh.vpos[current->next->vert]-mesh.vpos[current->vert]).norm2();
                            current = current->twin->next;
                        } while( current != half_0);

                        
                        k = 0;

                    }
                }
            }
        }
    }*/
#ifdef VERBOSE
    cout << endl << "simplified mesh from " << num_t_init << " to " << mesh.trigs.size() << " triangles." << endl;
#endif
}



// Here are some useful helper functions for the flip_edges function

int valence_cost(int valence) {
    int diff(valence - 6); // 6 is the optimal valence number
    return diff*diff;
}

vec3 cross_vec(vec3 const& a,vec3 const& b,vec3 const& c) {
    return (b-a).vec((c-b));
}

// Each edge of the mesh has two adjacent triangles. These two 
// triangles form together a shape with four corner vertices where
// two of them are connected to the edge. It is possible too, to 
// connect the edge to the other two vertices. If after such a change 
// the resulting edge is shorter, the overall valence number becomes 
// closer to 6, or the triangle shapes become more regular, it is useful 
// to change the connections of the edge accordingly. This is exactly what
// the following function does. the argument 'state' chooses between
// criterions for swapping the edge connections. There may be implemented
// further choices later on.
void flip_edges     (HalfedgeMesh& mesh, size_t state) {
    check_validity(mesh);
#ifdef VERBOSE
    cout << "FLIP-EDGES" << endl;
#endif
    /*
    // compute the valence of each vertex (number of neighbouring triangles/edges)
    vector<size_t> valences(mesh.verts.size(),0);

    for(Halfedge* elm : mesh.trigs) {
        valences[elm->vert]++;
        valences[elm->next->vert]++;
        valences[elm->next->next->vert]++;
    }
    

    size_t num_flipped(0);
    
    for(Halfedge* halfedge_A : mesh.edges){

        Halfedge* halfedge_B = halfedge_A->twin;

        if(halfedge_A == halfedge_B) {

        } else {

            size_t vert_0 = halfedge_A->vert;
            size_t vert_1 = halfedge_B->vert;

            size_t vert_A = halfedge_A->next->next->vert;
            size_t vert_B = halfedge_B->next->next->vert;
            

            // if valences are <=3, one must not reduce them further, otherwise
            // the mesh becomes invalid.
            if(valences[vert_0] > 3 and valences[vert_1] > 3) {

                int cost = 0;
                int cost_flip = 0;

                cost += valence_cost(valences[vert_0]);
                cost += valence_cost(valences[vert_1]);
                cost += valence_cost(valences[vert_A]);
                cost += valence_cost(valences[vert_B]);

                cost_flip += valence_cost(valences[vert_0]-1); // it is safe to subtract one here
                cost_flip += valence_cost(valences[vert_1]-1);
                cost_flip += valence_cost(valences[vert_A]+1);
                cost_flip += valence_cost(valences[vert_B]+1);


                // check new and old lengths:
                real length      = (mesh.vpos[vert_1] - mesh.vpos[vert_0]).norm2();
                real length_flip = (mesh.vpos[vert_A] - mesh.vpos[vert_B]).norm2();

                // It woule be possible too, to test the delaunay condition or other
                // shape quality proxys for triangles

                // test wether triangles are sufficiently similarly oriented
                vec3 nA = cross_vec(mesh.vpos[vert_0],mesh.vpos[vert_B],mesh.vpos[vert_A]);
                nA.normalize();
                vec3 nB = cross_vec(mesh.vpos[vert_B],mesh.vpos[vert_1],mesh.vpos[vert_A]);
                nB.normalize();
                real dotprod = nA.dot(nB);

                // condition whether or not to flip. Note that the variable state decides which criterion
                // will be applied
                bool flip = (state > 0 and cost_flip < cost) or (state == 0 and length_flip < length);

                // that the dotprod is larger than 0.8 is a necessary condition for flipping. The number 0.8 
                // was estimated through trials.
                if(flip and dotprod > 0.8) {
                    valences[vert_0]--;
                    valences[vert_1]--;
                    valences[vert_A]++;
                    valences[vert_B]++;

                    // make sure that the vertices 0 and 1 have valid pointers

                    mesh.verts[vert_0] = halfedge_B->next;
                    mesh.verts[vert_1] = halfedge_A->next;

                    // make sure that the triangles A and B have valid pointers

                    mesh.trigs[halfedge_A->trig] = halfedge_A;
                    mesh.trigs[halfedge_B->trig] = halfedge_B;

                    // flip the edge!

                    // changing the parameters of halfedge_A/B->next->next

                    Halfedge* temp_A = halfedge_A->next->next;
                    Halfedge* temp_B = halfedge_B->next->next;

                    temp_A->next = halfedge_B->next;
                    temp_B->next = halfedge_A->next;

                    // changing the parameters of halfedge_A/B->next

                    halfedge_A->next->next = halfedge_B;
                    halfedge_A->next->trig = halfedge_B->trig;

                    halfedge_B->next->next = halfedge_A;
                    halfedge_B->next->trig = halfedge_A->trig;

                    // changing the parameters of halfedge_A/B

                    halfedge_A->vert = temp_B->vert;
                    halfedge_B->vert = temp_A->vert;

                    halfedge_A->next = temp_A;
                    halfedge_B->next = temp_B;

                    num_flipped++;
                }
            }
        } 
    }*/
#ifdef VERBOSE
    cout << "flipped " << num_flipped << " edges." << endl;
#endif 
}

// this function puts each vertex to the average position
// of its direct neighbours. This has a smoothing character
void relax_vertices (Mesh& mesh) {
#ifdef VERBOSE
    cout << "RELAX-VERTICES" << endl;
#endif  
    /* 
    vector<vector<size_t>> neighbours = generate_neighbours(mesh);
    vector<vec3> new_vertices(mesh.verts.size());

    for(size_t i(0);i<mesh.verts.size();++i) {
        vec3 new_pos;
        real num(0.0);
        for(size_t j : neighbours[i]) {
            new_pos += mesh.verts[j];
            num++;
        }
        new_vertices[i] = (1.0/num)*new_pos;
    }

    mesh.verts = new_vertices;
    */
}

// the same function as above, but implemented for a HalfedgeMesh
void relax_vertices (HalfedgeMesh& mesh) {
#ifdef VERBOSE
    cout << "RELAX-VERTICES" << endl;
#endif
    /*
    vector<vec3> new_vertices(mesh.vpos.size());

    for(size_t i(0);i<mesh.verts.size();++i) {
        vec3 new_pos;
        real num(0.0);

        Halfedge* u(mesh.verts[i]);
        Halfedge* v = u;
        do {
            //if(v == v->twin) 
            v = v->twin;
            new_pos += mesh.vpos[v->vert];
            v = v->next;
            num++;
        } while(v != u);

        new_vertices[i] = (1.0/num)*new_pos;
    }

    mesh.vpos = new_vertices;
    */
}

// this function is used in project() and it finds the intersection of a ray, 
// starting from a position pos and pointing in direction dir, with a mesh.
// the result is given in the variables result and trig_index. If there are multiple 
// solutions, the one that is closest to pos is chosen. If there is no real solution,
// there will be still a solution returned which is not part of the Mesh. A corresponding
// test has to be implemented if needed. We suppose, that the position pos is close to the
// mesh surface and dir points into a similar direction as the normals in the vicinity of pos.
void trace_mesh(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index) {
    real s_min = -1.0;
    
    size_t m(mesh.trigs.size());

    for(size_t j(0);j<m;++j) {
        Triplet t(mesh.trigs[j]);

        vec3 a(mesh.verts[t.b]-mesh.verts[t.a]);
        vec3 b(mesh.verts[t.c]-mesh.verts[t.b]);
        vec3 c(mesh.verts[t.a]-mesh.verts[t.c]);

        vec3 n(a.vec(b));

        real s = n.dot(mesh.verts[t.a]-pos)/n.dot(dir);

        vec3 x = pos + s*dir;

        // find the position with minimal s - parameter
        if(abs(s) < s_min or s_min < 0.0) { 
            if(     (mesh.verts[t.a]-x).vec(a).dot(n) >= 0
                and (mesh.verts[t.b]-x).vec(b).dot(n) >= 0
                and (mesh.verts[t.c]-x).vec(c).dot(n) >= 0 ) {
                    s_min = abs(s);
                    result = x;
                    trig_index = j;
            }
        }
    }
}

// same function as above, but allowing only positive s. 
bool trace_mesh_positive(Mesh const& mesh,vec3 pos,vec3 dir,vec3& result,size_t& trig_index) {
    real s_min = 0.0;

    bool found(false);
    
    size_t m(mesh.trigs.size());

    for(size_t j(0);j<m;++j) {
        Triplet t(mesh.trigs[j]);

        vec3 a(mesh.verts[t.b]-mesh.verts[t.a]);
        vec3 b(mesh.verts[t.c]-mesh.verts[t.b]);
        vec3 c(mesh.verts[t.a]-mesh.verts[t.c]);

        vec3 n(a.vec(b));

        real s = n.dot(mesh.verts[t.a]-pos)/n.dot(dir);

        if(j == 0) s_min = s;

        vec3 x = pos + s*dir;

        // find the position with minimal s - parameter
        if(s > 0.0 and s < s_min) { 
            if(     (mesh.verts[t.a]-x).vec(a).dot(n) >= 0
                and (mesh.verts[t.b]-x).vec(b).dot(n) >= 0
                and (mesh.verts[t.c]-x).vec(c).dot(n) >= 0 ) {
                    s_min = s;
                    result = x;
                    trig_index = j;

                    found = true;
            }
        }
    }
    return found;
}

// this funciton projects all vertices of one mesh on the surface defined by
// the mesh 'other' by applying the function trace_mesh
void project(Mesh& mesh, Mesh const& other) {
    // normalized vertex normals
    vector<vec3> normals = generate_vertex_normals(mesh); // not nice!!

    size_t n(mesh.verts.size());

    vector<vec3> new_vertices(n);

    omp_set_num_threads(100);

    #pragma omp parallel
    {

    Mesh local_mesh,local_other;
    vector<vec3> local_normals;

    #pragma omp critical 
    {
        local_mesh = mesh;
        local_other = other;
        local_normals = normals;
    }


    #pragma omp for
    for(size_t i = 0;i<n;++i) {

        vec3 projected_pos;

        size_t index(0); // dummy variable

        trace_mesh(local_other,local_mesh.verts[i],local_normals[i],projected_pos,index);

        #pragma omp critical
        new_vertices[i] = projected_pos;
    }

    }

    mesh.verts.clear();
    mesh.verts = new_vertices;
}

// same function as above, but starting from origin and setting tracing-attempts without
// result to zero.
void project_from_origin(std::vector<vec3>& normals, Mesh const& other, real const& dist_to_wall) {
    
    size_t n(normals.size());

    omp_set_num_threads(100);

    #pragma omp parallel
    {

    Mesh local_other;
    vector<vec3> local_normals;

    #pragma omp critical 
    {
        local_normals = normals;
        local_other = other;
    }


    #pragma omp for
    for(size_t i = 0;i<n;++i) {

        vec3 projected_pos;

        size_t index(0); // dummmy variable

        bool found = trace_mesh_positive(local_other,vec3(),local_normals[i],projected_pos,index);

        #pragma omp critical 
        {
        
        if(found) {
            normals[i] = projected_pos
         } else {
            normals[i] = local_normals[i]*dist_to_wall;
            normals[i].x = 0.0;
         }

        }
    }

    }
}


void project_and_interpolate(Mesh& mesh, vector<real>& f_res, Mesh const& other, vector<real> const& f) {
    project_and_interpolate(mesh,generate_vertex_normals(mesh),f_res,other,f);
}

// with this function we can not only project one mesh on another, but also get the interpolated
// values of a piecewise linear function defined on the mesh 'other' by 'f' at the vertices of the 
// new mesh. Note too, that we are now not directly projecting on the flat triangles of the mesh 
// 'other', but on a local polynomial fit of the vertices of this mesh. This leads to a smoother 
// surface, especially in the case where the vertex density is locally much higher on the new mesh 
// than on the mesh 'other'.
void project_and_interpolate(Mesh& mesh,vector<vec3> const& vertex_normals, vector<real>& f_res, Mesh const& other, vector<real> const& f) {
#ifdef VERBOSE
    cout << "PROJECT-AND-INTERPOLATE" << endl;
#endif
    assert(f.size() == other.verts.size());
    vector<real> result(mesh.verts.size());

    // normalized vertex normals
    vector<vec3> normals = vertex_normals;

    vector<vec3> new_vertices(mesh.verts.size());

    // creating surface fits for each vertex of 'other'
    vector<vec3> other_normals = generate_vertex_normals(other);
    vector<vector<size_t>> verts = generate_neighbours(other); // alternatively generate_2_ring();
    vector<FittingTool> fits(other.verts.size());
    for(size_t i(0);i<other.verts.size();++i) {
        vector<vec3> positions;
        positions.push_back(other.verts[i]);
        for(size_t j : verts[i])
            positions.push_back(other.verts[j]);
        fits[i].compute_quadratic_fit(other_normals[i],other.verts[i],positions);
    }

    omp_set_num_threads(100);

    #pragma omp parallel
    {

    Mesh local_mesh,local_other;
    vector<vec3> local_normals;
    vector<vec3> local_other_normals;

    #pragma omp critical 
    {
        local_mesh = mesh;
        local_other = other;
        local_normals = normals;
        local_other_normals = other_normals;
    }

    size_t n(local_mesh.verts.size());

    #pragma omp for
    for(size_t i = 0;i<n;++i) {

        vec3 projected_pos;

        size_t index(0);

        trace_mesh(local_other,local_mesh.verts[i],local_normals[i],projected_pos,index);

        Triplet t = local_other.trigs[index];

        // copy coord system of the three fits at the three edges.
        // using transform, and then get_position to obtain the projection on the 
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

        // now we want to determine the coordinates of projected_pos
        // with respect to the usual triangle coordinates of the 
        // local triangle. These coordinates are then used for the 
        // interpolation of the scalar function and the weighting of
        // the fitting functions.

        // build local 2d coord system
        vec3 u(local_other.verts[t.b]-local_other.verts[t.a]);
        u.normalize();
        vec3 v(local_other.verts[t.c]-local_other.verts[t.b]);
        v = v - v.dot(u)*u;
        v.normalize();

        vec3 a(local_other.verts[t.b]-local_other.verts[t.a]);
        vec3 b(local_other.verts[t.c]-local_other.verts[t.b]);
        vec3 x(projected_pos         -local_other.verts[t.a]);

        real a0,a1,b0,b1,x0,x1,q,r;

        a0 = a.dot(u);
        a1 = a.dot(v);
        b0 = b.dot(u);
        b1 = b.dot(v);
        x0 = x.dot(u);
        x1 = x.dot(v);

        // q and r are the coordinates of x in local coordinates
        q = (x1*b0 - b1*x0)/(a1*b0-a0*b1); // q and r are well defined if a and b aren't collinear
        r = (a1*x0 - a0*x1)/(a1*b0-a0*b1);

        projected_pos = (1.0 - q)*Pa + (q - r)*Pb + r*Pc;

        #pragma omp critical 
        {
            new_vertices[i] = projected_pos;
            result[i] = (1.0 - q)*f[t.a] + (q - r)*f[t.b] + r*f[t.c];
        }
    }

    }

    mesh.verts.clear();
    mesh.verts = new_vertices;

    f_res.clear();
    f_res = result;
}

Mesh l2smooth(Mesh mesh,vector<size_t> const& vert_inds) {
    vector<vec3> normals = generate_vertex_normals(mesh);
    vector<vector<size_t>> ring = generate_2_ring(mesh);
    vector<vec3> copy = mesh.verts;
    for(size_t i : vert_inds) {
        vector<vec3> positions;
        positions.push_back(copy[i]);
        for(size_t j : ring[i])
            positions.push_back(copy[j]);

        FittingTool fit;
        fit.compute_quadratic_fit(normals[i],copy[i],positions);
        // refine the normal and recompute the fit (optional)
        vec3 ref_norm = fit.get_normal();
        fit.compute_quadratic_fit(ref_norm,copy[i],positions);

        mesh.verts[i] = fit.get_position(0.0,0.0);
    }

    return mesh;
}

Mesh l2smooth(Mesh mesh) {
    vector<size_t> inds(mesh.verts.size());
    iota(inds.begin(),inds.end(),0);
    return l2smooth(mesh,inds);
}

Mesh l2smooth(Mesh mesh, std::vector<real>& pot) {
    vector<size_t> inds(mesh.verts.size());
    iota(inds.begin(),inds.end(),0);
    return l2smooth(mesh,pot,inds);
}

Mesh l2smooth(Mesh mesh, std::vector<real>& pot, std::vector<size_t> const& vert_inds) {
    vector<vec3> normals = generate_vertex_normals(mesh);
    vector<vector<size_t>> ring = generate_2_ring(mesh);
    vector<vec3> copy = mesh.verts;
    vector<real> copy_pot = pot;
    for(size_t i : vert_inds) {
        vector<vec3> positions;
        vector<real> potential;
        positions.push_back(copy[i]);
        potential.push_back(copy_pot[i]);
        for(size_t j : ring[i]){
            positions.push_back(copy[j]);
            potential.push_back(copy_pot[j]);
        }

        FittingTool fit;
        fit.compute_quadratic_fit(normals[i],copy[i],positions);
        // refine the normal and recompute the fit (optional)
        vec3 ref_norm = fit.get_normal();
        fit.compute_quadratic_fit(ref_norm,copy[i],positions);

        mesh.verts[i] = fit.get_position(0.0,0.0);
        
        
        CoordSystem system(fit.copy_coord_system());
        for(size_t j(0);j<positions.size();++j) {
            vec3 trans = system.transform(positions[j]);
            trans.z = potential[j];
            positions[j] = trans;
        }

        fit.compute_quadratic_fit(vec3(0.0,0.0,1.0),vec3(),positions);
        vector<real> params(fit.get_params());

        pot[i] = params[0];
    }

    return mesh;
}

} // namespace Bem