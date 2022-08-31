
# Testing stuff in polygon.jl
let
    try
        p1 = ConvexPolygon([0,1,1,0], [0,0,1,1])
        # There should be an error!
        @test "Coordinates of polygons should be cyclical" == ""
    catch e
        @test 1==1
    end
    
    p1 = ConvexPolygon([0,1,1,0,0], [0,0,1,1,0])

    @test nvertices(p1) == 4
    @test area(p1) == 1
    @test normal(p1) > 0

    p2 = ConvexPolygon((0,0), (1,0), (1,1), (0,1), (0,0))

    @test nvertices(p2) == 4
    @test area(p2) == 1.0
    @test normal(p2) > 0
    @test centroid(p2) ≈ Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p2
    @test Point(0.5, -0.5) ∉ p2
    @test Point(eps(), eps()) ∈ p2
    @test Point(-eps(), -eps()) ∉ p2
    v2,conn2 = poly2mesh(p2)
    @test v2 == vertices(p2)
    @test conn2 == [1 2 3; 1 3 4]
    
    
    p3 = ConvexPolygon(reverse( [(0,0), (1,0), (1,1), (0,1), (0,0)] ))
    @test nvertices(p3) == 4
    @test area(p3) == 1.0
    @test normal(p3) < 0
    @test centroid(p3) == Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p3
    @test Point(0.5, -0.5) ∉ p3
    @test Point(eps(), eps()) ∈ p3
    @test Point(-eps(), -eps()) ∉ p3

    p4 = ConvexPolygon([(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,0)])
    @test nvertices(p4) == 4
    @test area(p4) == 1.0
    @test normal(p4) ≈ Vec(0.0, 0.0, 1.0)
    @test centroid(p4) ≈ Point(0.5, 0.5, 0.0)
    v4,conn4 = poly2mesh(p4)
    @test vertices(p4)==v4
    @test conn4 == [1 2 3; 1 3 4]
    
    
    p5 = ConvexPolygon(reverse([(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,0)]))
    
    @test area(p5) == 1.0
    @test normal(p5) == Vec(0.0, 0.0, -1.0)
    @test centroid(p4) == Point(0.5, 0.5, 0.0)
    
    
    
end
