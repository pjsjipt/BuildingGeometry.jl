using BuildingGeometry
using Test

#import GeometryBasics: Point, Point3, TriangleFace

@testset "BuildingGeometry.jl" begin

    include("test_polygon.jl")
    #include("test_polyhedron.jl")
    #include("test_polyhedronchop.jl")
    #include("test_discr_surface.jl")
    
end
