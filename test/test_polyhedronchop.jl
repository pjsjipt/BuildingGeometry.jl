# Testing functionality of polyhedronchop.jl

using LinearAlgebra

let B = BuildingGeometry

    p₀ = Point(0.0, 0.0, 0.0)
    n⃗ = Point(0.0, 0.0, 1.0)
    p₁ = Point(3.0, 1.0, 1.0)
    p₂ = Point(1.0, -1.0, -1.0)

    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ Point(2.0, 0.0, 0.0)

    p₀ = Point(1.0, 1.0, 1.0)
    n⃗ = Point(-1.0, 1.0, 0.0)
    p₁ = Point(0.5, 0.1, 0.1)
    p₂ = Point(0.5, 2.0, 0.1)
    
    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ Point(0.5, 0.5, 0.1)

    
    # Testing `cut_with_plane`

    p0 = Point(0.0, 0.0, 0.0)
    n = Point(-1.0, 0.0, 0.0)
    pts = Point3.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 4
    @test out[1] ≈ Point(0.0, -1.0, 0.0) 
    @test out[4] ≈ Point(0.0, 1.0, 0.0)
   
    pts = Point3.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
                   (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 6
    @test out[1] ≈ Point(0.0, -1.0, 0.0)
    @test out[3] ≈ Point(0.0, -0.5, 0.0)
    @test out[4] ≈ Point(0.0,  0.5, 0.0)
    @test out[6] ≈ Point(0.0,  1.0, 0.0)


    p0 = Point(0.0, 5.0, 5.0)
    n = Point(1.0, 0.0, 0.0)

    pts = Point3.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
                   (-1.0,1.0,0.0)])
    out = cut_with_plane(pts, p0, n, false)
    @test length(out) == 7
    @test out[1] == pts[1]
    @test out[2] ≈ Point(0.0, -1.0, 0.0)
    @test out[3] ≈ Point(0.0, -0.5, 0.0)
    @test out[4] ≈ Point(-1.0, 0.0, 0.0)
    @test out[5] ≈ Point(0.0,  0.5, 0.0)
    @test out[6] ≈ Point(0.0,  1.0, 0.0)
    @test out[7] ≈ Point(-1.0,  1.0, 0.0)

    out = cut_with_plane(pts, p0, n, true)
    @test length(out) == 7
    @test out[end] == pts[1]
    @test out[1] ≈ Point(0.0, -1.0, 0.0)
    @test out[2] ≈ Point(0.0, -0.5, 0.0)
    @test out[3] ≈ Point(-1.0, 0.0, 0.0)
    @test out[4] ≈ Point(0.0,  0.5, 0.0)
    @test out[5] ≈ Point(0.0,  1.0, 0.0)
    @test out[6] ≈ Point(-1.0,  1.0, 0.0)
    
    #@test out


    p0 = Point(0.0, 0.0, 0.0)
    n = Point(-1.0, 0.0, 0.0)
    pts = Point3.([(0.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0)])
    out = cut_with_plane(pts, p0, n, false)
    @test pts ≈ out

    out = cut_with_plane(pts, p0, n, true)
    @test pts[1] == out[end]

    
end


