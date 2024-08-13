
# Testing stuff in polygon.jl
using Unitful
import Unitful: m
m² = m*m
let
    
    p1 = ConvexPolygon([0,1,1,0], [0,0,1,1])

    @test nvertices(p1) == 4
    @test area(p1) == 1.0m²
    @test normal_(p1) > 0

    p2 = ConvexPolygon([(0,0), (1,0), (1,1), (0,1)])

    @test nvertices(p2) == 4
    @test area(p2) == 1.0m²
    @test normal_(p2) > 0
    @test centroid(p2) ≈ Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p2
    @test Point(0.5, -0.5) ∉ p2
    @test Point(eps(), eps()) ∈ p2
    @test Point(-eps(), -eps()) ∉ p2
    v2,conn2 = poly2mesh(p2)
    @test v2 == vertices(p2)
    @test conn2 == [1 2 3; 1 3 4]
    
    
    p3 = ConvexPolygon(reverse( [(0,0), (1,0), (1,1), (0,1)] ))
    @test nvertices(p3) == 4
    @test area(p3) == 1.0m²
    @test normal_(p3) < 0
    @test centroid(p3) == Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p3
    @test Point(0.5, -0.5) ∉ p3
    @test Point(eps(), eps()) ∈ p3
    @test Point(-eps(), -eps()) ∉ p3

    p4 = ConvexPolygon([(0,0,0), (1,0,0), (1,1,0), (0,1,0)])
    @test nvertices(p4) == 4
    @test area(p4) == 1.0m²
    @test normal(p4) ≈ Vec(0.0, 0.0, 1.0)
    @test centroid(p4) ≈ Point(0.5, 0.5, 0.0)
    v4,conn4 = poly2mesh(p4)
    @test vertices(p4)==v4
    @test conn4 == [1 2 3; 1 3 4]
    
    
    p5 = ConvexPolygon(reverse([(0,0,0), (1,0,0), (1,1,0), (0,1,0)]))
    
    @test area(p5) == 1.0m²
    @test normal(p5) == Vec(0.0, 0.0, -1.0)
    @test centroid(p4) == Point(0.5, 0.5, 0.0)
    

    p6 = ConvexPolygon([(0,0,0), (1,0,0), (0,4,0)])
    @test area(p6) ≈ 2.0m²
    @test centroid(p6) == Point(1/3, 4/3, 0.0)

    p7 = ConvexPolygon([(0,0,0), (1,0,0), (0,4,4)])
    @test area(p7) ≈ 2.0 * sqrt(2)m²
    @test centroid(p7) ≈ Point(1/3, 4/3, 4/3)
    @test normal(p7) ≈ Vec(0.0, -sqrt(2)/2, sqrt(2)/2)

    p8 = ConvexPolygon([(0,0,0), (1,0,0), (1,4,4), (0,4,4)])
    @test area(p8) ≈ 4.0 * sqrt(2)m²
    @test centroid(p8) ≈ Point(0.5, 2.0, 2.0)
    @test normal(p8) ≈ Vec(0.0, -sqrt(2)/2, sqrt(2)/2)
    
end
