module BuildingGeometry

using LinearAlgebra
import PyCall: pyimport, PyNULL 

const pyvoronoi = Ref(PyNULL())

function __init__()
    pyvoronoi[] =  pyimport("scipy.spatial").Voronoi
end



include("polygons.jl")
#include("polyhedron.jl")
#include("voronoi3d.jl")
#include("polyhedronchop.jl")
#include("discr_surface.jl")

end
