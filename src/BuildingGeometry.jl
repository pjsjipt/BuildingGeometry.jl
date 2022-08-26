module BuildingGeometry
import PyCall: pyimport, PyNULL 

const pyvoronoi = Ref(PyNULL())

function __init__()
    pyvoronoi[] =  pyimport("scipy.spatial").Voronoi
end



include("polygons.jl")
#include("surf_voronoi_project2d.jl")
include("polyhedron.jl")
include("voronoi3d.jl")
include("polyhedronchop.jl")

end
