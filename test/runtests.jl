using BuildingGeometry
using Test

import Meshes: Point, Vec, Point3, Triangle

@testset "BuildingGeometry.jl" begin

    include("test_polygon.jl")
    include("test_polyhedron.jl")
    include("test_polyhedronchop.jl")
    include("test_discr_surface.jl")
    
end
