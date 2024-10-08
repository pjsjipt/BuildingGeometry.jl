
# Testing stuff in polygon.jl
let
    
    p1 = ConvexPolygon([0.0,1,1,0], [0.0,0,1,1])

    @test nvertices(p1) == 4
    @test area(p1) == 1.0
    @test normalarea(p1) > 0

    p2 = ConvexPolygon([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)])

    @test nvertices(p2) == 4
    @test area(p2) == 1.0
    @test normalarea(p2) > 0
    @test centroid(p2) ≈ SVec(0.5, 0.5)
    @test SVec(0.5, 0.5) ∈ p2
    @test SVec(0.5, -0.5) ∉ p2
    @test SVec(eps(), eps()) ∈ p2
    @test SVec(-eps(), -eps()) ∉ p2
    v2,conn2 = poly2mesh(p2)
    @test v2 == vertices(p2)
    @test conn2 == [1 2 3; 1 3 4]
    
    
    p3 = ConvexPolygon(reverse( [(0.0,0.0), (1.0,0.0), (1.0,1), (0.0,1)] ))
    @test nvertices(p3) == 4
    @test area(p3) == 1.0
    @test normalarea(p3) < 0
    @test centroid(p3) == SVec(0.5, 0.5)
    @test SVec(0.5, 0.5) ∈ p3
    @test SVec(0.5, -0.5) ∉ p3
    @test SVec(eps(), eps()) ∈ p3
    @test SVec(-eps(), -eps()) ∉ p3

    p4 = ConvexPolygon([(0.0,0,0), (1.0,0,0), (1.0,1,0), (0.0,1,0)])
    @test nvertices(p4) == 4
    @test area(p4) == 1.0
    @test normal(p4) ≈ SVec(0.0, 0.0, 1.0)
    @test centroid(p4) ≈ SVec(0.5, 0.5, 0.0)
    v4,conn4 = poly2mesh(p4)
    @test vertices(p4)==v4
    @test conn4 == [1 2 3; 1 3 4]
    
    
    p5 = ConvexPolygon(reverse([(0.0,0,0), (1.0,0,0), (1.0,1,0), (0.0,1,0)]))
    
    @test area(p5) == 1.0
    @test normal(p5) == SVec(0.0, 0.0, -1.0)
    @test centroid(p4) == SVec(0.5, 0.5, 0.0)
    

    p6 = ConvexPolygon([(0.0,0,0), (1.0,0,0), (0.0,4,0)])
    @test area(p6) ≈ 2.0
    @test centroid(p6) == SVec(1/3, 4/3, 0.0)

    p7 = ConvexPolygon([(0.0,0,0), (1.0,0,0), (0.0,4,4)])
    @test area(p7) ≈ 2.0 * sqrt(2)
    @test centroid(p7) ≈ SVec(1/3, 4/3, 4/3)
    @test normal(p7) ≈ SVec(0.0, -sqrt(2)/2, sqrt(2)/2)

    p8 = ConvexPolygon([(0.0,0,0), (1.0,0,0), (1.0,4,4), (0.0,4,4)])
    @test area(p8) ≈ 4.0 * sqrt(2)
    @test centroid(p8) ≈ SVec(0.5, 2.0, 2.0)
    @test normal(p8) ≈ SVec(0.0, -sqrt(2)/2, sqrt(2)/2)
    
end
