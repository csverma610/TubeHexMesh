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
    Array3D  xyz;            // Coordinates of the node.
    MeshPtr  hexmesh;        // Hex mesh of the sphere.
    vector<EdgePtr> edges;   // edges attached to this node.
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

    std::array<NodePtr,2>  nodes;
    MeshPtr  hexmesh;
    std::vector<FacePtr>  profiles[2];   // Used for hexmesh...
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
        ofile << nodes.size() << " " << boundfaces.size() << " 0 " << endl;

        for( auto v : nodes)
            ofile << v->xyz[0] << " " << v->xyz[1] << " " << v->xyz[2] << endl;

        for( auto f : boundfaces)
            ofile << "4 " << f->nodes[0]->id << " "
                  << f->nodes[1]->id << " "
                  << f->nodes[2]->id << " "
                  << f->nodes[3]->id << endl;
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
                double x = -0.5 + i*dl;
                double y = -0.5 + j*dl;
                double z = -0.5 + k*dl;
                auto newnode = std::make_shared<Node>();
                newnode->xyz[0] = x;
                newnode->xyz[1] = y;
                newnode->xyz[2] = z;
                newnode->id     = index;
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
        v->xyz[0] = x;
        v->xyz[1] = y;
        v->xyz[2] = z;
        v->id     = i;
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
        newedge->nodes[0] = n0;
        newedge->nodes[1] = n1;
        graph->edges[i]   = newedge;
        n0->edges.push_back(newedge);
        n1->edges.push_back(newedge);
    }
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
	   newface->nodes[i] = newnode;
      }
      return newface;
}

//////////////////////////////////////////////////////////////////////////////

void buildTube( const NodePtr &vtx, const EdgePtr &edge)
{
    NodePtr othernode;

    int pos = 0;

    if( edge->nodes[0] == vtx) pos = 1;
    if( edge->nodes[1] == vtx) pos = 0;

    othernode = edge->nodes[pos];

    Array3D p0 = vtx->xyz;
    Array3D p1 = othernode->xyz;

    double dx  = p1[0] - p0[0];
    double dy  = p1[1] - p0[1];
    double dz  = p1[2] - p0[2];
    double dl  = sqrt(dx*dx + dy*dy + dz*dz);

    // A query point lies on the boundary of the junction.

    Array3D queryPoint;
    queryPoint[0] = vtx->radius*dx/dl;
    queryPoint[1] = vtx->radius*dy/dl;
    queryPoint[2] = vtx->radius*dz/dl;

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
    for( auto f : vtx->hexmesh->boundfaces)
    {
           double len = getLength(f->getCenter(), queryPoint);
	   if( len < mindist) {
               capface = f;
	       mindist   = len;
	   }
    }

    // Only one edge is supposed to pass through the face. If not, we
    // need to refine the junction node..
    assert(!capface->used_for_hex);

    capface->used_for_hex = 1;

    Array3D edgeCenter =  edge->getCenter();

    int nprofiles = 5;
    double dt = 1.0/(double)nprofiles;

    // Create profiles along the edge ( from the query point and edgecenter).

    Array3D p3d;
    edge->profiles[pos].push_back(capface);
    for( int i = 1; i < nprofiles-1; i++) {
         double t = i*dt;
	 p3d[0] = (1-t)*queryPoint[0] + t*edgeCenter[0];
	 p3d[1] = (1-t)*queryPoint[1] + t*edgeCenter[1];
	 p3d[2] = (1-t)*queryPoint[2] + t*edgeCenter[2];
	 auto newprofile = getProfileAt(capface, p3d);
	 edge->profiles[pos].push_back(newprofile);
    }

    // Create hex elements from the profile quads ....

    for( int i = 0; i < edge->profiles[pos].size()-1; i++) {
         auto face1 = edge->profiles[pos][i];
         auto face2 = edge->profiles[pos][i+1];
         auto cell = getHexElement(face1, face2);
	 edge->hexmesh->cells.push_back(cell);
    }
}
//////////////////////////////////////////////////////////////////////////////

void buildTube( const EdgePtr &edge)
{
    edge->hexmesh = std::make_shared<Mesh>();

    // Create tube starting from the node 0
    buildTube( edge->nodes[0], edge);

    // Create tube starting from the node 1
    buildTube( edge->nodes[1], edge);

    // Join the tube at the center of the edge. There is a danger of twisting
    // of the hex element. Need to address this issue.
    auto f1   = edge->profiles[0].back();
    auto f2   = edge->profiles[1].back();
    auto cell = getHexElement(f1,f2);
    edge->hexmesh->cells.push_back(cell);

    // Create hex mesh on the edge. Do not include the nodes of the junctions.
    // i.e. first profile on the either end. 

    for( int i = 1;  i < edge->profiles[0].size(); i++) {
	    auto f = edge->profiles[0][i];
	    for( int j  = 0; j < 4; j++) 
	    edge->hexmesh->nodes.push_back( f->nodes[j] );
    }

    for( int i = 1;  i < edge->profiles[1].size(); i++) {
	    auto f = edge->profiles[1][i];
	    for( int j  = 0; j < 4; j++) 
	    edge->hexmesh->nodes.push_back( f->nodes[j] );
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

    for( auto edge : vtx->edges) {
        auto n0 = edge->nodes[0];
        auto n1 = edge->nodes[1];
        double len = getLength( n0->xyz, n1->xyz);
        minlen  =  min(len, minlen);
    }
    assert( minlen > 1.0E-06);

    double radius = 0.01*minlen;

    vtx->hexmesh  = genSphere(radius);
    vtx->radius   = radius;
}

//////////////////////////////////////////////////////////////////////////////
MeshPtr accumulateHexMesh( const MeshPtr &graph)
{
   auto hexmesh = std::make_shared<Mesh>();

   for( auto v: graph->nodes) hexmesh->add(v->hexmesh);
   for( auto e: graph->edges) hexmesh->add(e->hexmesh);

   return hexmesh;
}
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    assert( argc == 2);

    auto  graph  = readGraph( argv[1] );

    cout << graph->edges.size() << endl;

    for( auto v: graph->nodes ) buildJunction(v);
    for( auto e: graph->edges ) {
	    buildTube(e);
    }

    auto hexmesh = accumulateHexMesh(graph);

    hexmesh->saveAs("hexmesh.off");

    return 0;
}
//////////////////////////////////////////////////////////////////////////////
