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


    pts = Point.([(1.0,-1.0,-1.0), (1.0,-1.0,+1.0), (1.0,1.0,1.0), (1.0,1.0,-1.0)])
    po = ConvexPolygon(pts)
    @test normal(po) ≈ Point(-4.0,0.0,0.0)
    n = Point(0.0, 0.0, 1.0)
    u = Point(1.0, 0.0, 0.0)
    v = Point(0.0, 1.0, 0.0)
    p₀ = Point(0.0, 0.0, 0.0)
    plane = B.Plane(p₀, n, u, v)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 2
    @test ipts[1] ≈ Point(1.0, -1.0, 0.0)
    @test ipts[2] ≈ Point(1.0, +1.0, 0.0)


    pts = Point.([(1.0,0.0,0.0), (1.0,+1.0,+1.0), (1.0,1.0,-1.0)])
    po = ConvexPolygon(pts)
    @test normal(po) ≈ Point(-1.0,0.0,0.0)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 2
    @test ipts[1] ≈ Point(1.0, 0.0, 0.0)
    @test ipts[2] ≈ Point(1.0, 1.0, 0.0)


    pts = Point.([(1.0,0.0,0.0), (1.0,+1.0,+1.0), (1.0,-1.0,1.0)])
    po = ConvexPolygon(pts)
    @test normal(po) ≈ Point(1.0,0.0,0.0)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 1
    @test ipts[1] ≈ Point(1.0, 0.0, 0.0)

    # Just circulate the points
    pts = Point.([(1.0,-1.0,1.0), (1.0,0.0,0.0), (1.0,+1.0,+1.0)])
    po = ConvexPolygon(pts)
    @test normal(po) ≈ Point(1.0,0.0,0.0)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 1
    @test ipts[1] ≈ Point(1.0, 0.0, 0.0)

    pts = Point.([(1.0,0.0, eps()), (1.0,+1.0,+1.0), (1.0,-1.0,1.0)])
    s = [(p-p₀)⋅n for p in pts]
    po = ConvexPolygon(pts)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 0

    pts = Point.([(1.0,0.0, -eps()), (1.0,+1.0,+1.0), (1.0,-1.0,1.0)])
    s = [(p-p₀)⋅n for p in pts]
    po = ConvexPolygon(pts)
    ipts = B.intersectfaceplane(po, plane)
    @test length(ipts) == 2
    
    
    
end


