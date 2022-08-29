using BuildingGeometry
using Test

import GeometryBasics: Point, Point3

@testset "BuildingGeometry.jl" begin

    include("test_polygon.jl")
    include("test_polyhedron.jl")
    include("test_polyhedronchop.jl")
    
end
