#pragma once

#include <string>
#include <array>
#include <vector>
#include <deque>
#include <set>
#include <fstream>

#include <QGLViewer/qglviewer.h>
#include <memory>

#include <QApplication>
#include <QKeyEvent>
#include <QColorDialog>

class Node;
typedef std::shared_ptr<Node> NodePtr;

class Edge;
typedef std::shared_ptr<Edge> EdgePtr;

class Face;
typedef std::shared_ptr<Face> FacePtr;

struct Node
{
    static NodePtr newObject();
    int  id;
    bool active   = 1;
    bool boundary = 0;
    bool visit    = 0;
    std::array<float,3> xyz;
    std::vector<EdgePtr> edges;
    std::vector<FacePtr> faces;
};

inline double length2( const NodePtr n0, const NodePtr &n1)
{
    double dx = n0->xyz[0] - n1->xyz[0];
    double dy = n0->xyz[1] - n1->xyz[1];
    double dz = n0->xyz[2] - n1->xyz[2];

    double l2 = dx*dx + dy*dy + dz*dz;
    return l2;
}

inline NodePtr Node:: newObject()
{
    NodePtr v(new Node);
    return v;
}

struct Edge {
    static EdgePtr newObject( NodePtr &n0, NodePtr &n1);

    Edge() {};
    Edge( const NodePtr &v0, const NodePtr &v1) {
        nodes[0] = v0;
        nodes[1] = v1;
    }

    bool hasNodes( const NodePtr &n0, const NodePtr &n1) const {
        if( (nodes[0]->id == n0->id) && (nodes[1]->id == n1->id) ) return 1;
        if( (nodes[0]->id == n1->id) && (nodes[1]->id == n0->id) ) return 1;
        return 0;
    }

    bool isBoundary() const {
        if( faces[1] == nullptr) return 1;
    }

    bool active    = 1;
    bool interface = 0;
    std::array<NodePtr,2>  nodes    = {nullptr,nullptr};
    std::array<FacePtr,2>  faces = {nullptr,nullptr};
};

inline EdgePtr Edge:: newObject(NodePtr &n0, NodePtr &n1)
{
    EdgePtr e(new Edge(n0,n1));
    return e;
}

struct Face {
    static FacePtr newObject( NodePtr &n0, NodePtr &n1, NodePtr &n2, NodePtr &n3);


    Face() {};
    Face( NodePtr &v0, NodePtr &v1, NodePtr &v2, NodePtr &v3) {
        nodes[0] = v0;
        nodes[1] = v1;
        nodes[2] = v2;
        nodes[3] = v3;
    }

    std::array<float,3> getCentroid() const;
    float getArea() const;

    bool active  = 1;
    bool visited = 0;
    int  id;
    std::array<NodePtr,4> nodes;
    std::array<EdgePtr,4> edges;
};

inline std::array<float,3> Face :: getCentroid() const
{
    std::array<float,3> center = {0.0, 0.0, 0.0};

    for( int i = 0; i < 4; i++) {
        auto p = nodes[i]->xyz;
        center[0] += p[0];
        center[1] += p[1];
        center[2] += p[2];
    }
    center[0] /= 3.0;
    center[1] /= 3.0;
    center[2] /= 3.0;

    return center;
}

inline FacePtr Face:: newObject(NodePtr &n0, NodePtr &n1, NodePtr &n2, NodePtr &n3)
{
    auto f = std::make_shared<Face>(n0,n1,n2,n3);
    return f;
}

struct Mesh
{
    EdgePtr addEdge( NodePtr &n0, NodePtr &n1, FacePtr &f);
    void    addFace( FacePtr &f);

    std::vector<NodePtr> nodes;
    std::vector<EdgePtr> edges;
    std::vector<FacePtr> faces;

    double radius;
    std::array<double,3> center = {0.0, 0.0, 0.0};


    void saveAs( const std::string &s);
};


class QuadViewer: public QGLViewer
{
public:

    void readMesh( const std::string &s);
    void readAffinityMatrix( const std::string &s);

protected:
    virtual void draw();
    virtual void init();
    virtual void keyPressEvent( QKeyEvent *e);
    virtual void mousePressEvent( QMouseEvent *e);
    virtual void mouseReleaseEvent( QMouseEvent *e);

private:
    Mesh srcmesh, dstmesh, currmesh;
    std::array<double,3> startPos, endPos;
    int  maxSteps = 100;
    std::array<float,3>  bgColor = {0.2, 0.2, 0.2};

    int  pickEntity = 0;
    bool displayWires = 1;
    bool displaySurface = 1;
    bool useLights      = 0;
    bool displayIDs     = 0;

    void drawFaces(const Mesh &themesh);
    void drawEdges(const Mesh &themesh);
};
