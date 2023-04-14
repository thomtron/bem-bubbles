#include "MeshIO.hpp"

#include <iostream>
#include <fstream> // for ofstream
#include <algorithm> // for clamp
#include <array>

using namespace std;


namespace Bem {

void export_obj (string filename, Mesh const& mesh)
{
    ofstream output(filename);

    output << "# output from Mesh class" << endl;
    
    for(vec3 const& v : mesh.verts){
        output << "v " << v.x << ' ' << v.y << ' ' << v.z << endl;
    }

    // note in .obj the indices start with 1 - we thus have to add 1 to each index.
    for(Triplet const& t : mesh.trigs){
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

template<typename T>
struct Vertex5 {
    T x,y,z,v,w;
    Vertex5(T x,T y,T z,T v,T w)
        :x(x),y(y),z(z),v(v),w(w) {}
};

void export_ply (string filename, Mesh const& mesh) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

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
        float_vec vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z);
        output.write((char*) (&vec),sizeof(float_vec));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply (string filename, Mesh const& mesh, vector<real> const& values, real min, real max) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

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
        float_vec vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z);
        output.write((char*) (&vec),sizeof(float_vec));
        char red(clamp((values[i]-min)/(max-min)*255,0.,255.));
        char green(0);
        char blue(clamp((max-values[i])/(max-min)*255,0.,255.));
        output.write(&red,sizeof(char));
        output.write(&green,sizeof(char));
        output.write(&blue,sizeof(char));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply_colors (string filename, Mesh const& mesh, vector<vec3> const& colors) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

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
        float_vec vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z);
        output.write((char*) (&vec),sizeof(float_vec));
        char red    = clamp(colors[i].x*255.,0.,255.);
        char green  = clamp(colors[i].y*255.,0.,255.);
        char blue   = clamp(colors[i].z*255.,0.,255.);
        output.write(&red,sizeof(char));
        output.write(&green,sizeof(char));
        output.write(&blue,sizeof(char));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}


void export_ply_float (string filename, Mesh const& mesh, vector<real> const& values) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

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
        Vertex<float> vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z,values[i]);
        output.write((char*) (&vec),sizeof(Vertex<float>));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply_double (string filename, Mesh const& mesh, vector<real> const& values) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

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
        Vertex<double> vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z,values[i]);
        output.write((char*) (&vec),sizeof(Vertex<double>));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply_float (string filename, Mesh const& mesh, vector<real> const& phi, vector<real> const& psi) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property float x" << endl; 
    output << "property float y" << endl; 
    output << "property float z" << endl;
    output << "property float phi" << endl;
    output << "property float psi" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        Vertex5<float> vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z,phi[i],psi[i]);
        output.write((char*) (&vec),sizeof(Vertex5<float>));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply_double (string filename, Mesh const& mesh, vector<real> const& phi, vector<real> const& psi) {
    ofstream output(filename);

    size_t n(mesh.verts.size());
    size_t m(mesh.trigs.size());

    output << "ply" << endl;
    output << "format binary_little_endian 1.0" << endl;
    output << "comment Output from BEM simulation" << endl;
    output << "element vertex " << n << endl;
    output << "property double x" << endl; 
    output << "property double y" << endl; 
    output << "property double z" << endl;
    output << "property double phi" << endl;
    output << "property double psi" << endl;
    output << "element face " << m << endl;
    output << "property list uchar uint vertex_indices" << endl;
    output << "end_header" << endl;

    output.close();

    output.open(filename, std::ofstream::binary | std::ofstream::app);

    for(size_t i(0);i<n;++i) {
        Vertex5<double> vec(mesh.verts[i].x,mesh.verts[i].y,mesh.verts[i].z,phi[i],psi[i]);
        output.write((char*) (&vec),sizeof(Vertex5<double>));
    }
    for(Triplet const& t : mesh.trigs){
        uint_vec vec(t.a,t.b,t.c);
        char num(3);
        output.write(&num,sizeof(char));
        output.write((char*) (&vec),sizeof(uint_vec));
    }

    output.close();
}

void export_ply_float_separat (string filename, Mesh const& mesh, vector<real> const& values) {
    ofstream output(filename);

    size_t m(mesh.trigs.size());

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

    output.open(filename, ofstream::binary | ofstream::app);

    for(size_t i(0);i<m;++i) {
        for(size_t j(0);j<3;++j){
            size_t k = mesh.trigs[i][j];
            Vertex<float> vec(mesh.verts[k].x,mesh.verts[k].y,mesh.verts[k].z,values[3*i+j]);
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

void export_ply_double_separat (string filename, Mesh const& mesh, vector<real> const& values) {
    ofstream output(filename);

    size_t m(mesh.trigs.size());

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
            size_t k = mesh.trigs[i][j];
            Vertex<double> vec(mesh.verts[k].x,mesh.verts[k].y,mesh.verts[k].z,values[3*i+j]);
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



struct float_5 {
    float x,y,z,v,w;
};

struct double_5 {
    double x,y,z,v,w;
};

struct float_4 {
    float x,y,z,v;
};

struct double_4 {
    double x,y,z,v;
};

struct float_3 {
    float x,y,z;
};

struct double_3 {
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

array<size_t,7> propSize {
    sizeof(float),
    sizeof(double),
    sizeof(unsigned int),
    sizeof(int),
    sizeof(char),
    sizeof(unsigned int),
    sizeof(int)
};


vector<Property> read_properties(ifstream& input) {
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


void import_ply (std::string filename, Mesh& mesh, std::vector<real>& values) {

    // first, clean everything up.
    mesh.verts.clear();
    mesh.trigs.clear();
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
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                    values.push_back(result.v);
                }
            }
            else {
#ifdef VERBOSE
                cout << "reading 3 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_3 result = parse_binary<float_3>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                }
            }
        }
        else if(vertexprops[0] == DOUBLE and vertexprops[1] == DOUBLE and vertexprops[2] == DOUBLE and vertexprops[3] == DOUBLE){
#ifdef VERBOSE
            cout << "reading 4 doubles" << endl;
#endif
            for(size_t i(0);i<num_vertex;++i){
                double_4 result = parse_binary<double_4>(binarystream,vertex_size);
                mesh.verts.push_back(vec3(result.x,result.y,result.z));
                values.push_back(result.v);
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
                    mesh.trigs.push_back(Triplet(result.i,result.j,result.k));
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
                    mesh.trigs.push_back(Triplet(result.i,result.j,result.k));
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



void import_ply (string filename, Mesh& mesh, vector<real>& phi, vector<real>& psi) {

    // first, clean everything up.
    mesh.verts.clear();
    mesh.trigs.clear();
    phi.clear();
    psi.clear();


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
            if(vertexprops[3] == FLOAT and vertexprops[4] == FLOAT){
#ifdef VERBOSE
                cout << "reading 5 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_5 result = parse_binary<float_5>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                    phi.push_back(result.v);
                    psi.push_back(result.w);
                }
            }
            else if(vertexprops[3] == FLOAT) {
#ifdef VERBOSE
                cout << "reading 4 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_4 result = parse_binary<float_4>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                    phi.push_back(result.v);
                }
            }
            else {
#ifdef VERBOSE
                cout << "reading 3 floats" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    float_3 result = parse_binary<float_3>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                }
            }
        }
        else if(vertexprops[0] == DOUBLE and vertexprops[1] == DOUBLE and vertexprops[2] == DOUBLE){
            if(vertexprops[3] == DOUBLE and vertexprops[4] == DOUBLE){
#ifdef VERBOSE
                cout << "reading 5 doubles" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    double_5 result = parse_binary<double_5>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                    phi.push_back(result.v);
                    psi.push_back(result.w);
                }
            }
            else if(vertexprops[3] == DOUBLE) {
#ifdef VERBOSE
                cout << "reading 4 doubles" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    double_4 result = parse_binary<double_4>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                    phi.push_back(result.v);
                }
            }
            else {
#ifdef VERBOSE
                cout << "reading 3 doubles" << endl;
#endif
                for(size_t i(0);i<num_vertex;++i){
                    double_3 result = parse_binary<double_3>(binarystream,vertex_size);
                    mesh.verts.push_back(vec3(result.x,result.y,result.z));
                }
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
                    mesh.trigs.push_back(Triplet(result.i,result.j,result.k));
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
                    mesh.trigs.push_back(Triplet(result.i,result.j,result.k));
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


void import_ply(std::string filename, Mesh& mesh) {
    vector<Bem::real> vals;
    import_ply(filename,mesh,vals);
}



} // namespace Bem