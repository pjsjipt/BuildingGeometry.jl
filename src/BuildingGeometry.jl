module BuildingGeometry

using LinearAlgebra
import PyCall: pyimport, PyNULL 

const pyvoronoi = Ref(PyNULL())

function __init__()
    pyvoronoi[] =  pyimport("scipy.spatial").Voronoi
end

import Meshes
import Meshes: Manifold, 𝔼, CRS, ucross, udot, unormalize
import Meshes: Point, area, centroid, normal, coords, Ring, vertices, Polygon
import Meshes: nvertices, Vec, isclosed, measure, Triangle, Box
import Meshes: Polyhedron, measure, volume, Point, nvertices, vertices, nfacets
import Meshes: boundingbox, Plane, to, connect, SimpleMesh, discretize, simplexify
using Unitful
import Unitful: ustrip, unit

export Point, normalarea, normal_
export ConvexPolygon, area, centroid, normal, coordinates, vertices, nvertices
export poly2mesh, measure, volume, nfacets, boundingbox
export ConvexPolyhedron
export cut_with_plane, chopwithpolyhedron, slicemesh, intersectmesh!, intersect_tri
export NodeInfo, nodearea, nodepoint, nodeside, nodetag
export BuildingSurface, savebsurf, loadbsurf, buildsurface, buildingslice, mergemeshes
export readraw,tri2mesh
export addforcecontrib!, forcematrix
export reescalemesh, translatemesh, rotatemesh


include("polygons.jl")
include("polyhedron.jl")
include("voronoi3d.jl")
#include("chopmesh.jl")
#include("discr_surface.jl")
#include("node_info.jl")
#include("intersect_mesh.jl")
#include("building_surface.jl")
#include("raw.jl")
#include("forces.jl")
#include("transform.jl")

end
