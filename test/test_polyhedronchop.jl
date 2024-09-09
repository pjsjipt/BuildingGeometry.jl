# Testing functionality of polyhedronchop.jl

using LinearAlgebra

let B = BuildingGeometry

    p₀ = SVec(0.0, 0.0, 0.0)
    n⃗ = SVec(0.0, 0.0, 1.0)
    p₁ = SVec(3.0, 1.0, 1.0)
    p₂ = SVec(1.0, -1.0, -1.0)

    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ SVec(2.0, 0.0, 0.0)

    p₀ = SVec(1.0, 1.0, 1.0)
    n⃗ = SVec(-1.0, 1.0, 0.0)
    p₁ = SVec(0.5, 0.1, 0.1)
    p₂ = SVec(0.5, 2.0, 0.1)
    
    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ SVec(0.5, 0.5, 0.1)

    
    # Testing `cut_with_plane`

    p0 = SVec(0.0, 0.0, 0.0)
    n = SVec(-1.0, 0.0, 0.0)
    pts = SVec.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 4
    @test out[1] ≈ SVec(0.0, -1.0, 0.0) 
    @test out[4] ≈ SVec(0.0, 1.0, 0.0)
   
    pts = SVec.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
                   (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 6
    @test out[1] ≈ SVec(0.0, -1.0, 0.0)
    @test out[3] ≈ SVec(0.0, -0.5, 0.0)
    @test out[4] ≈ SVec(0.0,  0.5, 0.0)
    @test out[6] ≈ SVec(0.0,  1.0, 0.0)


    p0 = SVec(0.0, 5.0, 5.0)
    n = SVec(1.0, 0.0, 0.0)

    pts = SVec.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
                  (-1.0,1.0,0.0)])
    out = cut_with_plane(pts, p0, n, false)
    @test length(out) == 7
    @test out[1] == pts[1]
    @test out[2] ≈ SVec(0.0, -1.0, 0.0)
    @test out[3] ≈ SVec(0.0, -0.5, 0.0)
    @test out[4] ≈ SVec(-1.0, 0.0, 0.0)
    @test out[5] ≈ SVec(0.0,  0.5, 0.0)
    @test out[6] ≈ SVec(0.0,  1.0, 0.0)
    @test out[7] ≈ SVec(-1.0,  1.0, 0.0)

    out = cut_with_plane(pts, p0, n, true)
    @test length(out) == 7
    @test out[end] == pts[1]
    @test out[1] ≈ SVec(0.0, -1.0, 0.0)
    @test out[2] ≈ SVec(0.0, -0.5, 0.0)
    @test out[3] ≈ SVec(-1.0, 0.0, 0.0)
    @test out[4] ≈ SVec(0.0,  0.5, 0.0)
    @test out[5] ≈ SVec(0.0,  1.0, 0.0)
    @test out[6] ≈ SVec(-1.0,  1.0, 0.0)

    p0 = SVec(0.0, 0.0, 0.0)
    n = SVec(-1.0, 0.0, 0.0)
    pts = SVec.([(0.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0)])
    out = cut_with_plane(pts, p0, n, false)
    @test all(pts .≈ out)

    out = cut_with_plane(pts, p0, n, true)
    @test pts[1] == out[end]


    # Testing the chopping itsel

    pts0 = makebbox( (-1.0, 1.0), (-2.0, 2.0), (-4.0, 4.0) )
    pts  = SVec.(pts0...)

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)


    # Trivial case: no intersection, outside the bounding boxes
    tri = Tri(SVec(10.0,0.0,0.0), SVec(10.5,0.0,0.0), SVec(0.0,10.5,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 0

    # Trivial case:  intersection, triangle inside polyhedron
    tri = Tri(SVec(0.0,0.0,0.0), SVec(0.5,0.0,0.0), SVec(0.0,0.5,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 1
    @test vertices(ch[1]) == vertices(tri)

    tri = Tri(SVec(-2.0,0.0,0.0), SVec(2.0,-2.0,0.0), SVec(2.0,2.0,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 2
    # The exact nodes are tricky.
    pts = [SVec(-1.0, -0.5, 0.0), SVec(1.0, -1.5, 0.0), SVec(1.0, 1.5, 0.0),
           SVec(-1.0, 0.5, 0.0)]
    # I know the solution is two triangles
    ii = 0
    for (k,pt) in enumerate(pts)
        if pt==vertices(ch[1])[1]
            ii = k
        end
    end
#    @test vertices(ch[1])[2] ≈ pts[4]
#    @test vertices(ch[1])[3] ≈ pts[1]

#    @test vertices(ch[2])[2] ≈ pts[1]
#    @test vertices(ch[2])[3] ≈ pts[2]
    
end


