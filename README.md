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

