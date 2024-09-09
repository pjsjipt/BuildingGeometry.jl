
# Testing polyhedron.jl stuff

let
    
    pts0 = makebbox( (0.0, 1.0), (0.0, 2.0), (0.0, 4.0) )
    pts  = SVec.(pts0...)

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)
    pf = [ConvexPolygon(pts[ff]) for ff in faces]

    @test nfacets(pp) == 6
    @test nvertices(pp) == 8
    
    for i in 1:length(faces)
        @test pp[i].contour == pf[i].contour
    end
    
    @test area(pp[1]) == 4.0
    @test area(pp[2]) == 8.0
    @test area(pp[3]) == 4.0
    @test area(pp[4]) == 8.0
    @test area(pp[5]) == 2.0
    @test area(pp[6]) == 2.0

    @test normal(pp[1]) ≈ SVec( 0.0,-1.0, 0.0)
    @test normal(pp[2]) ≈ SVec( 1.0, 0.0, 0.0)
    @test normal(pp[3]) ≈ SVec( 0.0, 1.0, 0.0)
    @test normal(pp[4]) ≈ SVec(-1.0, 0.0, 0.0)
    @test normal(pp[5]) ≈ SVec( 0.0, 0.0,-1.0) 
    @test normal(pp[6]) ≈ SVec( 0.0, 0.0, 1.0)
    
    @test volume(pp) == 8.0
    @test centroid(pp) ≈ SVec(0.5, 1.0, 2.0)

    @test SVec(0.5, 0.5, 0.5) ∈ pp
    @test SVec(0.5, 0.5,-0.5) ∉ pp
    @test SVec(eps(), eps(), eps()) ∈ pp 
    @test SVec(-eps(), eps(), eps()) ∉ pp



    pts  = SVec.([(0.0,0,0), (1.0,0.0,0), (0.0,2,0), (1.0,2,0),
                   (4.0,0,4), (5.0,0,4), (4.0,2,4), (5.0,2,4)])

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)
    pf = [ConvexPolygon(pts[ff]) for ff in faces]

    
    @test area(pp[1]) == 4.0
    @test area(pp[2]) == 8.0 * sqrt(2)
    @test area(pp[3]) == 4.0
    @test area(pp[4]) == 8.0 * sqrt(2)
    @test area(pp[5]) == 2.0
    @test area(pp[6]) == 2.0
    s2 = sqrt(2)/2
    @test normal(pp[1]) ≈ SVec( 0.0,-1.0, 0.0)
    @test normal(pp[2]) ≈ SVec( s2, 0.0, -s2)
    @test normal(pp[3]) ≈ SVec( 0.0, 1.0, 0.0)
    @test normal(pp[4]) ≈ SVec(-s2, 0.0, s2)
    @test normal(pp[5]) ≈ SVec( 0.0, 0.0,-1.0) 
    @test normal(pp[6]) ≈ SVec( 0.0, 0.0, 1.0)
    
    @test volume(pp) == 8.0
    @test centroid(pp) ≈ SVec(2.5, 1.0, 2.0)

    @test SVec(0.5, 0.5, 0.25) ∈ pp
    @test SVec(0.5, 0.5,-0.5) ∉ pp
    @test SVec(2*eps(), eps(), eps()) ∈ pp 
    @test SVec(eps(), eps(), 2*eps()) ∉ pp 
    @test SVec(0.5, 1.0, 0.5 + eps()) ∉ pp
    @test SVec(0.5, 1.0, 0.5 - eps()) ∈ pp
    
    
end
