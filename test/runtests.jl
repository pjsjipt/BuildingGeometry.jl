using BuildingGeometry
using Test

import Meshes: Point, Vec,  Triangle

@testset "BuildingGeometry.jl" begin

    include("test_polygon.jl")
    include("test_polyhedron.jl")
    include("test_polyhedronchop.jl")
    include("test_discr_surface.jl")
    #include("test_intersect_mesh.jl")
    #include("test_buildsurface.jl")
    #include("test_forces.jl")
end
