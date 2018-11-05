#include "QuadViewer.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
void QuadViewer::init()
{
    setSelectRegionWidth(10);
    setSelectRegionHeight(10);
}
////////////////////////////////////////////////////////////////////////////////

EdgePtr Mesh:: addEdge( NodePtr &n0, NodePtr &n1, FacePtr &face)
{
    NodePtr vmin = std::min(n0,n1);

    for( auto oldedge : vmin->edges) {
        if( oldedge->hasNodes(n0,n1)) {
            oldedge->faces[1] = face;
            return oldedge;
        }
    }

    EdgePtr newedge(new Edge(n0,n1));
    newedge->faces[0] = face;
    vmin->edges.push_back(newedge);
    edges.push_back(newedge);
    return newedge;
}
////////////////////////////////////////////////////////////////////////////////

void Mesh::addFace( FacePtr &newface)
{
    assert( newface );
    auto n0 = newface->nodes[0]; assert( n0 );
    auto n1 = newface->nodes[1]; assert( n1 );
    auto n2 = newface->nodes[2]; assert( n2 );
    auto n3 = newface->nodes[3]; assert( n3 );
    assert((n0 != n1) && (n1 != n2) && (n2 != n0));

    newface->edges[0] = addEdge(n0,n1,newface);
    newface->edges[1] = addEdge(n1,n2,newface);
    newface->edges[2] = addEdge(n2,n3,newface);
    newface->edges[3] = addEdge(n3,n0,newface);

    n0->faces.push_back(newface);
    n1->faces.push_back(newface);
    n2->faces.push_back(newface);
    n3->faces.push_back(newface);

    faces.push_back(newface);
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer:: readMesh( const string &filename)
{
    ifstream ifile( filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Warning: Input file not read " << endl;
        return;
    }
    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Warning: Input file not in Off format" << endl;
        return;
    }

    size_t numNodes, numFaces, numEdges;
    ifile >> numNodes >> numFaces >> numEdges;

    double x, y, z;

    srcmesh.nodes.resize(numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        ifile >> x >> y >> z;
        NodePtr v = Node::newObject();
        v->xyz[0] = x;
        v->xyz[1] = y;
        v->xyz[2] = z;
        v->id     = i;
        srcmesh.nodes[i] = v;
    }

    size_t index = 0;
    int dummy, v0, v1, v2, v3;

    for( size_t i = 0; i < numFaces; i++) {
        ifile >> dummy >> v0 >> v1 >> v2 >> v3;
        assert( dummy == 4);
        NodePtr n0   = srcmesh.nodes[v0];
        NodePtr n1   = srcmesh.nodes[v1];
        NodePtr n2   = srcmesh.nodes[v2];
        NodePtr n3   = srcmesh.nodes[v3];
        FacePtr newface = Face::newObject(n0,n1,n2,n3);
        newface->id     = i;
        srcmesh.addFace(newface);
    }
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer::keyPressEvent( QKeyEvent *e)
{
	/*
    if( e->key() == Qt::Key_0) {
        pickEntity = 0;
        this->setSelectedName(-1);
    }

    if( e->key() == Qt::Key_N) {
        displayIDs = !displayIDs;
    }

    if( e->key() == Qt::Key_S) {
        displaySurface = !displaySurface;
        update();
        return;
    }

    if( e->key() == Qt::Key_W) {
        displayWires = !displayWires;
    }

    if( e->key() == Qt::Key_L) {
        useLights = !useLights;
        update();
        return;
    }

    if( e->key() == Qt::Key_N) {
        nstep++;
        double t = nstep*dt;
        if( t <= 1.0) {
            At  = AffineLib::expSE(t*logA);
            mult( At, currmesh);
            update();
        }
        return;
    }

    if( e->key() == Qt::Key_R) {
	    nstep = 1;
            double t = nstep*dt;
            if( t <= 1.0) {
            At  = AffineLib::expSE(t*logA);
            mult( At, currmesh);
            update();
            }
            return;
    }
    if( e->key() == Qt::Key_Home) {
        qglviewer::Vec pos;
        pos[0]  = srcmesh.center[0];
        pos[1]  = srcmesh.center[1];
        pos[2]  = srcmesh.center[2];
        camera()->setSceneCenter(pos);
        camera()->setSceneRadius(srcmesh.radius);
        camera()->centerScene();
        camera()->showEntireScene();
        update();
        return;
    }
    */

    QGLViewer::keyPressEvent(e);


    update();
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer:: mousePressEvent( QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer::mouseReleaseEvent( QMouseEvent *e)
{
    int id = this->selectedName();

    QGLViewer::mouseReleaseEvent(e);

    update();
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer::drawFaces(Mesh &themesh)
{
    if( useLights ) glEnable(GL_LIGHTING);

    for ( auto f : themesh.faces) {
        if(f->active) {
            glBegin(GL_QUADS);
            glVertex3fv( &f->nodes[0]->xyz[0] );
            glVertex3fv( &f->nodes[1]->xyz[0] );
            glVertex3fv( &f->nodes[2]->xyz[0] );
            glVertex3fv( &f->nodes[3]->xyz[0] );
            glEnd();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void QuadViewer::draw()
{
    glPolygonOffset(1.0,1.0);
    glEnable(GL_POLYGON_OFFSET_LINE);

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);
    glColor3f( 1.0, 0.0, 0.0);
    drawFaces(srcmesh);

    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
    glColor3f( 0.0, 0.0, 0.0);
    drawFaces(srcmesh);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  assert( argc == 2);
  QApplication application(argc, argv);

  // Instantiate the viewer.
  QuadViewer viewer;

  viewer.setWindowTitle("QuadViewer");
  viewer.readMesh( argv[1] );

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
}

