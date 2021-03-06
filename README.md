# TubeHexMesh

In general, automatic hex mesh generation is a hard problem. But for some simple
geometries, it is possible to generate high quality hex mesh. Here we provide
a completely automatic hexmesh for a graph given by the nodes and edges.

Following the "TetWild" philosophy, we first generate a toplogically correct mesh
and geometry later. With this design, we generate a perfect hexmesh at each 
node and connect edges by tubes. An edge is connected by two nodes and from each
node a tube in generated and these two tubes are glued at the center of the
edge.

There may be many challenges in improving the geometric qualities of the mesh
(1) Elements may be twisted (2) Global intersections (3) Concave eells (4)
degenerate elements (i.e. volume zero) etc.  Further research is required to improve
the hexmesh quality.


## Example: xross.off

OFF  
5 0 4 
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
-1.0 0.0 0.0
0.0 -1.0 0.0
0 1
0 2
0 3
0 4

Should produce Hexmesh as:

![alt text](./example.png "Title")


## If problems occcur

1.   If a node has many edges then try refining the node hemesh. i.e. Change the second parameter "n" in genSphere(radius, n)
     However, remember that with the size "n", 6*(n-1)*(n-1) edges can be supported. So increment n only when it is 
     absolutely necessary.
2.   A tube is twisted at the center. This will occur when two profiles are incorrectly combined to form a hex element.
     Try changing the ordering of the nodes of one of the profiles.



