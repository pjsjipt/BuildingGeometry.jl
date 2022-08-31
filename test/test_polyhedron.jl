
# Testing polyhedron.jl stuff

let

    pts0 = makebbox( (0.0, 1.0), (0.0, 2.0), (0.0, 4.0) )
    pts  = Point.(pts0...)

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)
    pf = [ConvexPolygon(pts[[ff;ff[1]]]) for ff in faces]

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

    @test normal(pp[1]) ≈ Vec( 0.0,-1.0, 0.0)
    @test normal(pp[2]) ≈ Vec( 1.0, 0.0, 0.0)
    @test normal(pp[3]) ≈ Vec( 0.0, 1.0, 0.0)
    @test normal(pp[4]) ≈ Vec(-1.0, 0.0, 0.0)
    @test normal(pp[5]) ≈ Vec( 0.0, 0.0,-1.0) 
    @test normal(pp[6]) ≈ Vec( 0.0, 0.0, 1.0)
    
    @test volume(pp) == 8.0
    @test centroid(pp) ≈ Point(0.5, 1.0, 2.0)

    @test Point(0.5, 0.5, 0.5) ∈ pp
    @test Point(0.5, 0.5,-0.5) ∉ pp
    @test Point(eps(), eps(), eps()) ∈ pp 
    @test Point(-eps(), eps(), eps()) ∉ pp



    pts  = Point.([(0,0,0), (1,0,0), (0,2,0), (1,2,0),
                   (4,0,4), (5,0,4), (4,2,4), (5,2,4)])

    faces = [[1,2,6,5], [2,4,8,6], [4,3,7,8], [1,5,7,3], [1,3,4,2], [5,6,8,7]]
    pp = ConvexPolyhedron(pts, faces)
    pf = [ConvexPolygon(pts[[ff;ff[1]]]) for ff in faces]

    
    @test area(pp[1]) == 4.0
    @test area(pp[2]) == 8.0 * sqrt(2)
    @test area(pp[3]) == 4.0
    @test area(pp[4]) == 8.0 * sqrt(2)
    @test area(pp[5]) == 2.0
    @test area(pp[6]) == 2.0
    s2 = sqrt(2)/2
    @test normal(pp[1]) ≈ Vec( 0.0,-1.0, 0.0)
    @test normal(pp[2]) ≈ Vec( s2, 0.0, -s2)
    @test normal(pp[3]) ≈ Vec( 0.0, 1.0, 0.0)
    @test normal(pp[4]) ≈ Vec(-s2, 0.0, s2)
    @test normal(pp[5]) ≈ Vec( 0.0, 0.0,-1.0) 
    @test normal(pp[6]) ≈ Vec( 0.0, 0.0, 1.0)
    
    @test volume(pp) == 8.0
    @test centroid(pp) ≈ Point(2.5, 1.0, 2.0)

    @test Point(0.5, 0.5, 0.25) ∈ pp
    @test Point(0.5, 0.5,-0.5) ∉ pp
    @test Point(2*eps(), eps(), eps()) ∈ pp 
    @test Point(eps(), eps(), 2*eps()) ∉ pp 
    @test Point(0.5, 1.0, 0.5 + eps()) ∉ pp
    @test Point(0.5, 1.0, 0.5 - eps()) ∈ pp
    
    
end
