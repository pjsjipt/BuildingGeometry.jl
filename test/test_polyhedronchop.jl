# Testing functionality of polyhedronchop.jl

using LinearAlgebra

let B = BuildingGeometry

    p₀ = Point(0.0, 0.0, 0.0)
    n⃗ = Vec(0.0, 0.0, 1.0)
    p₁ = Point(3.0, 1.0, 1.0)
    p₂ = Point(1.0, -1.0, -1.0)

    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ Point(2.0, 0.0, 0.0)

    p₀ = Point(1.0, 1.0, 1.0)
    n⃗ = Vec(-1.0, 1.0, 0.0)
    p₁ = Point(0.5, 0.1, 0.1)
    p₂ = Point(0.5, 2.0, 0.1)
    
    p = B.intersectpoint(n⃗, p₀, p₁, p₂)
    @test p ≈ Point(0.5, 0.5, 0.1)

    
    # Testing `cut_with_plane`

    p0 = Point(0.0, 0.0, 0.0)
    n = Vec(-1.0, 0.0, 0.0)
    pts = Point.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 4
    @test out[1] ≈ Point(0.0, -1.0, 0.0) 
    @test out[4] ≈ Point(0.0, 1.0, 0.0)
   
    pts = Point.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
                   (-1.0,1.0,0.0)])
    
    out = cut_with_plane(pts, p0, n)
    @test length(out) == 6
    @test out[1] ≈ Point(0.0, -1.0, 0.0)
    @test out[3] ≈ Point(0.0, -0.5, 0.0)
    @test out[4] ≈ Point(0.0,  0.5, 0.0)
    @test out[6] ≈ Point(0.0,  1.0, 0.0)


    p0 = Point(0.0, 5.0, 5.0)
    n = Vec(1.0, 0.0, 0.0)

    pts = Point.([(-1.0,-1.0,0.0), (1.0,-1.0,0.0), (-1.0, 0.0, 0.0), (1.0,1.0,0.0),
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

    p0 = Point(0.0, 0.0, 0.0)
    n = Vec(-1.0, 0.0, 0.0)
    pts = Point3.([(0.0,-1.0,0.0), (1.0,-1.0,0.0), (1.0,1.0,0.0), (0.0,1.0,0.0)])
    out = cut_with_plane(pts, p0, n, false)
    @test all(pts .≈ out)

    out = cut_with_plane(pts, p0, n, true)
    @test pts[1] == out[end]


    # Testing the chopping itsel

    pts0 = makebbox( (-1.0, 1.0), (-2.0, 2.0), (-4.0, 4.0) )
    pts  = Point.(pts0...)

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)


    # Trivial case: no intersection, outside the bounding boxes
    tri = TriangleFace(Point(10.0,0.0,0.0), Point(10.5,0.0,0.0), Point(0.0,10.5,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 0

    # Trivial case:  intersection, triangle inside polyhedron
    tri = TriangleFace(Point(0.0,0.0,0.0), Point(0.5,0.0,0.0), Point(0.0,0.5,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 1
    @test coordinates(ch[1]) == coordinates(tri)

    tri = TriangleFace(Point(-2.0,0.0,0.0), Point(2.0,-2.0,0.0), Point(2.0,2.0,0.0))
    ch  = chopwithpolyhedron(pp, tri)
    @test length(ch) == 2
    # The exact nodes are tricky.
    pts = [Point(-1.0, -0.5, 0.0), Point(1.0, -1.5, 0.0), Point(1.0, 1.5, 0.0),
           Point(-1.0, 0.5, 0.0)]
    # I know the solution is two triangles

    @test area(pts) ≈ sum(area.(ch))
    
    p1 = ConvexPolygon(pts)
    numeq = [0,0]
    for (i,t) in enumerate(ch)
        for v in coordinates(t)
            for p in pts
                if v ≈ p
                    numeq[i] += 1
                end
            end
        end
    end
    @test numeq[1] == 3
    @test numeq[2] == 3
            
    
    #@test coordinates(ch[1])[2] ≈ pts[4]
    #@test coordinates(ch[1])[3] ≈ pts[1]

    #@test coordinates(ch[2])[2] ≈ pts[1]
    #@test coordinates(ch[2])[3] ≈ pts[2]
    
end


