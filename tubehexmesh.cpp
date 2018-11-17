#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>
#include <limits>
#include <map>
#include <set>

using namespace std;

typedef std::array<double,3> Array3D;
typedef std::array<int,4>    Array4I;
typedef std::array<int,8>    Array8I;

class Node;
typedef std::shared_ptr<Node>  NodePtr;

class Edge;
typedef std::shared_ptr<Edge>  EdgePtr;

class Face;
typedef std::shared_ptr<Face>  FacePtr;

class Mesh;
typedef std::shared_ptr<Mesh>  MeshPtr;

struct Node
{
    size_t   id;             // Unique ID of the node.
    double   radius;         // Radius of the sphere.
    bool     fixed  = 0;     // Can it move or not?
    Array3D  xyz;            // Coordinates of the node.
    MeshPtr  hexmesh;        // Hex mesh of the sphere.
    vector<EdgePtr> beams;   // beam attached to this node.
    vector<EdgePtr> meshedges;
};

struct Edge
{
    Array3D getCenter() const {
        Array3D pc = {0.0, 0.0, 0.0};
        for( int i = 0; i < 2; i++) {
            pc[0] += nodes[i]->xyz[0];
            pc[1] += nodes[i]->xyz[1];
            pc[2] += nodes[i]->xyz[2];
        }
        pc[0] *= 0.50;
        pc[1] *= 0.50;
        pc[2] *= 0.50;
        return pc;
    }

    double getLength() const
    {
        double dx = nodes[1]->xyz[0] - nodes[0]->xyz[0];
        double dy = nodes[1]->xyz[1] - nodes[0]->xyz[1];
        double dz = nodes[1]->xyz[2] - nodes[0]->xyz[2];
        return sqrt(dx*dx + dy*dy + dz*dz);
    }

    std::array<NodePtr,2>  nodes;
    MeshPtr  hexmesh;
    std::vector<FacePtr>  profiles[2];   // Used for hexmesh...
    std::vector<FacePtr>  quadprofiles;
    std::array<FacePtr,2> capfaces;
};

///////////////////////////////////////////////////////////////////////

struct Face
{
    std::array<NodePtr,4>  nodes;  // Quad-face
    Array3D getCenter() const {
        Array3D pc = {0.0, 0.0, 0.0};
        for( int i = 0; i < 4; i++) {
            pc[0] += nodes[i]->xyz[0];
            pc[1] += nodes[i]->xyz[1];
            pc[2] += nodes[i]->xyz[2];
        }
        pc[0] *= 0.25;
        pc[1] *= 0.25;
        pc[2] *= 0.25;
        return pc;
    }

    bool used_for_hex = 0;

};

///////////////////////////////////////////////////////////////////////

struct Cell
{
    std::array<NodePtr,8>  nodes;

    FacePtr getFaceAt(int dir)
    {
        auto face = std::make_shared<Face>();
        switch( dir)
        {
        case -1:
            face->nodes = { nodes[0], nodes[4], nodes[7], nodes[3] };
            break;
        case +1:
            face->nodes = { nodes[1], nodes[2], nodes[6], nodes[5] };
            break;
        case -2:
            face->nodes = { nodes[0], nodes[1], nodes[5], nodes[4] };
            break;
        case +2:
            face->nodes = { nodes[3], nodes[7], nodes[6], nodes[2] };
            break;
        case -3:
            face->nodes = { nodes[0], nodes[3], nodes[2], nodes[1] };
            break;
        case +3:
            face->nodes = { nodes[4], nodes[5], nodes[6], nodes[7] };
            break;
        }
        return face;
    }
};

typedef std::shared_ptr<Cell>  CellPtr;
double sphRadius;

///////////////////////////////////////////////////////////////////////

struct Mesh
{
    std::vector<NodePtr> nodes;
    std::vector<EdgePtr> edges;
    std::vector<CellPtr> cells;

    std::vector<NodePtr> boundnodes;
    std::vector<FacePtr> boundfaces;

    void add( const MeshPtr &m) {
        for( auto v: m->nodes) nodes.push_back(v);
        for( auto c: m->cells) cells.push_back(c);
    }

    void saveAs( const string &filename) {

        ofstream ofile( filename.c_str(), ios::out);

        ofile << "OFF" << endl;
        ofile << nodes.size() << " " << cells.size() << " 0 " << endl;

        for( auto v : nodes)
            ofile << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;

        for( auto c : cells)
            ofile << "8 " << c->nodes[0]->id << " "
                  << c->nodes[1]->id << " "
                  << c->nodes[2]->id << " "
                  << c->nodes[3]->id << " "
                  << c->nodes[4]->id << " "
                  << c->nodes[5]->id << " "
                  << c->nodes[6]->id << " "
                  << c->nodes[7]->id << endl;
    }

    void print( ofstream &ofile, FacePtr &f)
    {
        ofile << "4 " << f->nodes[0]->id << " "
              << f->nodes[1]->id << " "
              << f->nodes[2]->id << " "
              << f->nodes[3]->id << endl;
    };

    void saveFaces( const string &filename) {

        ofstream ofile( filename.c_str(), ios::out);

        ofile << "OFF" << endl;
        ofile << nodes.size() << " " << 6*cells.size() << " 0 " << endl;

        for( auto v : nodes)
            ofile << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;

        for( auto c : cells) {
            auto f0 = c->getFaceAt(-1);
            print( ofile, f0);
            auto f1 = c->getFaceAt(+1);
            print( ofile, f1);
            auto f2 = c->getFaceAt(-2);
            print( ofile, f2);
            auto f3 = c->getFaceAt(+2);
            print( ofile, f3);
            auto f4 = c->getFaceAt(-3);
            print( ofile, f4);
            auto f5 = c->getFaceAt(+3);
            print( ofile, f5);
        }
    }
};

///////////////////////////////////////////////////////////////////////

CellPtr getHexElement( const FacePtr &face1, const FacePtr &face2)
{
    auto newcell = std::make_shared<Cell>();

    newcell->nodes[0] = face1->nodes[0];
    newcell->nodes[1] = face1->nodes[1];
    newcell->nodes[2] = face1->nodes[2];
    newcell->nodes[3] = face1->nodes[3];

    newcell->nodes[4] = face2->nodes[0];
    newcell->nodes[5] = face2->nodes[1];
    newcell->nodes[6] = face2->nodes[2];
    newcell->nodes[7] = face2->nodes[3];

    return newcell;
}

///////////////////////////////////////////////////////////////////////

MeshPtr genSphere(double radius = 1.0, int n = 4)
{
    double dl = 2*radius/(double)(n-1);

    auto newmesh = std::make_shared<Mesh>();
    newmesh->nodes.resize(n*n*n);

    auto isBoundary = []( int i, int j, int k, int n) {
        if( i == 0 || i == n-1) return 1;
        if( j == 0 || j == n-1) return 1;
        if( k == 0 || k == n-1) return 1;
        return 0;
    };

    auto getOffset = []( int i, int j, int k, int n) {
        return k*n*n + j*n + i;
    };

    int index = 0;
    for( int k = 0; k < n; k++) {
        for( int j = 0; j < n; j++) {
            for( int i = 0; i < n; i++) {
                double x = -radius + i*dl;
                double y = -radius + j*dl;
                double z = -radius + k*dl;
                auto newnode = std::make_shared<Node>();
                newnode->xyz = {x,y,z};
                newnode->id  = index;
                newmesh->nodes[index++] = newnode;
                if( isBoundary(i,j,k,n) ) newmesh->boundnodes.push_back(newnode);
            }
        }
    }

    newmesh->cells.resize( (n-1)*(n-1)*(n-1));
    index = 0;
    int  offset;
    std::array<NodePtr,8> cellnodes;
    for( int k = 0; k < n-1; k++) {
        for( int j = 0; j < n-1; j++) {
            for( int i = 0; i < n-1; i++) {
                offset = k*n*n + j*n + i;
                cellnodes[0] = newmesh->nodes[offset];
                cellnodes[1] = newmesh->nodes[offset+1];
                cellnodes[2] = newmesh->nodes[offset+n+1];
                cellnodes[3] = newmesh->nodes[offset+n];

                offset += n*n;
                cellnodes[4] = newmesh->nodes[offset];
                cellnodes[5] = newmesh->nodes[offset+1];
                cellnodes[6] = newmesh->nodes[offset+n+1];
                cellnodes[7] = newmesh->nodes[offset+n];

                auto newcell = std::make_shared<Cell>();
                newcell->nodes = cellnodes;
                newmesh->cells[index++] = newcell;
                if( i == 0   ) newmesh->boundfaces.push_back( newcell->getFaceAt(-1));
                if( i == n-2 ) newmesh->boundfaces.push_back( newcell->getFaceAt(+1));
                if( j == 0   ) newmesh->boundfaces.push_back( newcell->getFaceAt(-2));
                if( j == n-2 ) newmesh->boundfaces.push_back( newcell->getFaceAt(+2));
                if( k == 0   ) newmesh->boundfaces.push_back( newcell->getFaceAt(-3));
                if( k == n-2 ) newmesh->boundfaces.push_back( newcell->getFaceAt(+3));
            }
        }
    }

    double normal;
    for( auto v : newmesh->boundnodes) {
        auto xyz  = v->xyz;
        double dl = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
        double nx = xyz[0]/dl;
        double ny = xyz[1]/dl;
        double nz = xyz[2]/dl;
        v->xyz[0] = radius*nx;
        v->xyz[1] = radius*ny;
        v->xyz[2] = radius*nz;
    }

    return newmesh;
}
//////////////////////////////////////////////////////////////////////////////

void translate( MeshPtr &mesh, const Array3D &p3d)
{
    for( auto v : mesh->nodes) {
        v->xyz[0] += p3d[0];
        v->xyz[1] += p3d[1];
        v->xyz[2] += p3d[2];
    }
}

//////////////////////////////////////////////////////////////////////////////

MeshPtr readGraph( const string &filename)
{
    ifstream ifile( filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Warning: Input file not read " << endl;
        return nullptr;
    }
    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Warning: Input file not in Off format" << endl;
        return nullptr;
    }

    size_t numNodes, numFaces, numEdges;
    ifile >> numNodes >> numFaces >> numEdges;

    assert( numFaces == 0);

    double x, y, z;

    auto graph = std::make_shared<Mesh>();

    graph->nodes.resize(numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        ifile >> x >> y >> z;
        auto v = std::make_shared<Node>();
        v->xyz = {x,y,z};
        v->id  = i;
        graph->nodes[i] = v;
    }

    size_t index = 0;
    int dummy, v0, v1;

    graph->edges.resize(numEdges);
    for( size_t i = 0; i < numEdges; i++) {
        ifile >> v0 >> v1;
        auto n0 = graph->nodes[v0];
        auto n1 = graph->nodes[v1];
        auto newedge = std::make_shared<Edge>();
        newedge->nodes  = {n0,n1};
        graph->edges[i] = newedge;
        n0->beams.push_back(newedge);
        n1->beams.push_back(newedge);
    }


    auto getLength = [] ( const Array3D &p0, const Array3D &p1)
    {
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];

        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    double minlen = std::numeric_limits<double>::max();

    for( auto e : graph->edges) {
        auto n0 = e->nodes[0];
        auto n1 = e->nodes[1];
        double len = getLength( n0->xyz, n1->xyz);
        minlen  =  min(len, minlen);
    }
    assert( minlen > 1.0E-06);

    sphRadius = 0.05*minlen;

    return graph;
}

//////////////////////////////////////////////////////////////////////////////

FacePtr getProfileAt(const FacePtr &refface, const Array3D &pos)
{
    // Make a copy of the refface and translate to given position...

    Array3D refpos = refface->getCenter();
    double  dx  = pos[0] - refpos[0];
    double  dy  = pos[1] - refpos[1];
    double  dz  = pos[2] - refpos[2];

    auto newface = std::make_shared<Face>();

    for( int i = 0; i < 4; i++) {
        auto newnode = std::make_shared<Node>();
        newnode->xyz[0] = refface->nodes[i]->xyz[0] + dx;
        newnode->xyz[1] = refface->nodes[i]->xyz[1] + dy;
        newnode->xyz[2] = refface->nodes[i]->xyz[2] + dz;
        newnode->fixed  = 1;
        newface->nodes[i] = newnode;
    }
    return newface;
}

//////////////////////////////////////////////////////////////////////////////
void translate( const FacePtr &face, const Array3D &newCenter)
{
    auto currCenter = face->getCenter();
    double dx = newCenter[0] - currCenter[0];
    double dy = newCenter[1] - currCenter[1];
    double dz = newCenter[2] - currCenter[2];

    for( int i = 0; i < 4; i++) {
        face->nodes[i]->xyz[0] += dx;
        face->nodes[i]->xyz[1] += dy;
        face->nodes[i]->xyz[2] += dz;
    }
}
//////////////////////////////////////////////////////////////////////////////

void findCap( const EdgePtr &edge, int pos)
{
    auto v0 = edge->nodes[pos];
    auto v1 = edge->nodes[(pos+1)%2];

    auto p0 = v0->xyz;
    auto p1 = v1->xyz;

    double dx  = p1[0] - p0[0];
    double dy  = p1[1] - p0[1];
    double dz  = p1[2] - p0[2];
    double dl  = sqrt(dx*dx + dy*dy + dz*dz);

    // A query point lies on the boundary of the junction.

    Array3D queryPoint;
    queryPoint[0] = v0->radius*dx/dl + v0->xyz[0];
    queryPoint[1] = v0->radius*dy/dl + v0->xyz[1];
    queryPoint[2] = v0->radius*dz/dl + v0->xyz[2];

    auto getLength = [] ( const Array3D &p0, const Array3D &p1)
    {
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];

        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    double mindist = std::numeric_limits<double>::max();

    // Check will boundary face intersects the query point.

    FacePtr capface;
    for( auto f : v0->hexmesh->boundfaces)
    {
        double len = getLength(f->getCenter(), queryPoint);
        if( len < mindist) {
            capface = f;
            mindist = len;
        }
    }

    // Only one edge is supposed to pass through the face. If not, we
    // need to refine the junction node..
    assert(!capface->used_for_hex);

    capface->used_for_hex = 1;
    edge->capfaces[pos] = capface;

//  translate( capface, queryPoint);
}

//////////////////////////////////////////////////////////////////////////////
void match( FacePtr &face, const FacePtr &refface)
{
    auto getLength = [] ( const Array3D &p0, const Array3D &p1)
    {
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];

        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    int pos0 = -1;

    double minlen = std::numeric_limits<double>::max();

    // look for the closest node on face to the "0th" node on refface ..

    Array3D p0, p1, p2;
    p0 = refface->nodes[0]->xyz;
    for( int j = 0; j < 4; j++) {
        p1  = face->nodes[j]->xyz;
        auto len = getLength(p0,p1);
        if( len < minlen) {
            pos0 = j;
            minlen = len;
        }
    }

    // Now look for the closest node on face to the 1st node on refface.
    // This node can not the node already found and it also cannot the third
    // node.

    minlen = std::numeric_limits<double>::max();

    int pos1 = -1;
    p0 = refface->nodes[1]->xyz;
    for( int j = 0; j < 4; j++) {
        if( (j != pos0) && (j != (pos0+2)%4) ) {
            p1  = face->nodes[j]->xyz;
            auto len = getLength(p0,p1);
            if( len < minlen) {
                pos1 = j;
                minlen = len;
            }
        }
    }

    // Now reorder the face with
    std::array<NodePtr,4> newnodes;
    newnodes[0] = face->nodes[pos0];
    newnodes[1] = face->nodes[pos1];
    newnodes[2] = face->nodes[(pos0+2)%4];
    newnodes[3] = face->nodes[(pos1+2)%4];

    face->nodes = newnodes;
}

////////////////////////////////////////////////////////////////////////////////

void buildTube( const EdgePtr &edge)
{
    findCap( edge, 0);
    findCap( edge, 1);

    int nprofiles =  max(5.0, edge->getLength()/sphRadius);

    double dt = 1.0/(double)(nprofiles-1);

    edge->quadprofiles.clear();

    auto face0 = edge->capfaces[0];
    auto face1 = edge->capfaces[1];

    assert( face0 != face1);

    edge->quadprofiles.push_back(face0);

    edge->hexmesh = std::make_shared<Mesh>();

    Array3D xyz;
    for( int i = 2; i < nprofiles-1; i++) {
        double t = i*dt;
        xyz[0] = (1-t)*edge->nodes[0]->xyz[0] + t*edge->nodes[1]->xyz[0];
        xyz[1] = (1-t)*edge->nodes[0]->xyz[1] + t*edge->nodes[1]->xyz[1];
        xyz[2] = (1-t)*edge->nodes[0]->xyz[2] + t*edge->nodes[1]->xyz[2];
        auto newface = getProfileAt( face0, xyz);
        edge->hexmesh->nodes.push_back(newface->nodes[0] );
        edge->hexmesh->nodes.push_back(newface->nodes[1] );
        edge->hexmesh->nodes.push_back(newface->nodes[2] );
        edge->hexmesh->nodes.push_back(newface->nodes[3] );
        edge->quadprofiles.push_back(newface);
    }

//   match(face1, edge->quadprofiles.back());

    edge->quadprofiles.push_back( face1 );

    for( int i = 0; i < edge->quadprofiles.size()-1; i++) {
        auto face1 = edge->quadprofiles[i];
        auto face2 = edge->quadprofiles[i+1];
        auto cell  = getHexElement(face1, face2);
        edge->hexmesh->cells.push_back(cell);
    }
}

//////////////////////////////////////////////////////////////////////////////

void buildJunction( const NodePtr &vtx)
{
    auto getLength = [] ( const Array3D &p0, const Array3D &p1)
    {
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];

        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    double minlen = std::numeric_limits<double>::max();

    for( auto edge : vtx->beams) {
        auto n0 = edge->nodes[0];
        auto n1 = edge->nodes[1];
        double len = getLength( n0->xyz, n1->xyz);
        minlen  =  min(len, minlen);
    }
    assert( minlen > 1.0E-06);

    double radius = 0.01*minlen;

    vtx->hexmesh  = genSphere(sphRadius,5);
    vtx->radius   = sphRadius;
    translate( vtx->hexmesh, vtx->xyz);
}
//
//////////////////////////////////////////////////////////////////////////////
//
void smoothJunction( const NodePtr &vtx)
{
    std::map<NodePtr, set<NodePtr>> vmap;
    for( auto f: vtx->hexmesh->boundfaces) {
        for( int i = 0; i < 4; i++) {
            auto n0 = f->nodes[i];
            auto n1 = f->nodes[(i+1)%4];
            vmap[n0].insert(n1);
            vmap[n1].insert(n0);
        }
    }

    Array3D xyz;

    int pos = 0;
    FacePtr f0, f1;
    for( auto e : vtx->beams) {
        if( e->nodes[0] == vtx ) pos = 0;
        if( e->nodes[1] == vtx ) pos = 1;
        int nprofiles = e->quadprofiles.size();
        if( pos == 0) {
            f0 = e->quadprofiles[0];
            f1 = e->quadprofiles[1];
        } else {
            f0 = e->quadprofiles[nprofiles-1];
            f1 = e->quadprofiles[nprofiles-2];
        }
        for( int i = 0; i < 4; i++) {
            auto n0 = f0->nodes[i];
            auto n1 = f1->nodes[i];
            vmap[n0].insert(n1);
            vmap[n1].insert(n0);
        }
    }

    std::map<NodePtr, Array3D> newPos;

    for( int i = 0; i < 5; i++) {
        for( auto keyval : vmap) {
            auto vi = keyval.first;
            if( !vi->fixed) {
                auto vset = keyval.second;
                Array3D xyz = {0.0, 0.0, 0.0};
                for( auto vj : vset ) {
                    xyz[0] += vj->xyz[0];
                    xyz[1] += vj->xyz[1];
                    xyz[2] += vj->xyz[2];
                }
                newPos[vi][0] = xyz[0]/(double)vmap[vi].size();
                newPos[vi][1] = xyz[1]/(double)vmap[vi].size();
                newPos[vi][2] = xyz[2]/(double)vmap[vi].size();
            }
        }

        for( auto keyval : newPos ) {
            auto vi   = keyval.first;
            vi->xyz   = keyval.second;
        }

        double normal;
        auto center = vtx->xyz;
        for( auto v: vtx->hexmesh->boundnodes) {
            xyz[0]  = v->xyz[0] - center[0];
            xyz[1]  = v->xyz[1] - center[1];
            xyz[2]  = v->xyz[2] - center[2];
            double dl = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
            double nx = xyz[0]/dl;
            double ny = xyz[1]/dl;
            double nz = xyz[2]/dl;
            v->xyz[0] = sphRadius*nx + center[0];
            v->xyz[1] = sphRadius*ny + center[1];
            v->xyz[2] = sphRadius*nz + center[2];
        }
    }
}
//
//////////////////////////////////////////////////////////////////////////////
//
MeshPtr accumulateHexMesh( const MeshPtr &graph)
{
    auto hexmesh = std::make_shared<Mesh>();

    for( auto v: graph->nodes) hexmesh->add(v->hexmesh);
    for( auto e: graph->edges) hexmesh->add(e->hexmesh);

    int index = 0;
    for( auto v : hexmesh->nodes) v->id = index++;

    return hexmesh;
}
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    assert( argc == 2);

    auto  graph  = readGraph( argv[1] );

    for( auto v: graph->nodes ) buildJunction(v);
    for( auto e: graph->edges ) buildTube(e);

    for( auto v: graph->nodes ) smoothJunction(v);

    auto hexmesh = accumulateHexMesh(graph);

    hexmesh->saveAs("hexmesh.off");
    hexmesh->saveFaces("facemesh.off");

    return 0;
}
//////////////////////////////////////////////////////////////////////////////
