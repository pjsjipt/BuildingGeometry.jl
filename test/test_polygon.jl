
# Testing stuff in polygon.jl
import GeometryBasics: Point
let
    p1 = ConvexPolygon([0,1,1,0], [0,0,1,1])

    @test length(coordinates(p1)) == 4
    @test area(p1) == 1
    @test normal(p1) > 0

    p2 = ConvexPolygon(Point{2,Float64}.([(0,0), (1,0), (1,1), (0,1)]))

    @test length(coordinates(p2)) == 4
    @test area(p2) == 1.0
    @test normal(p2) > 0
    @test centroid(p2) == Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p2
    @test Point(0.5, -0.5) ∉ p2
    @test Point(eps(), eps()) ∈ p2
    @test Point(-eps(), -eps()) ∉ p2
    v2,conn2 = poly2mesh(p2)
    @test v2 == coordinates(p2)
    @test conn2 == [1 2 3; 1 3 4]
    
    
    p3 = ConvexPolygon(Point{2,Float64}.( reverse( [(0,0), (1,0), (1,1), (0,1)] )))
    @test length(coordinates(p3)) == 4
    @test area(p3) == 1.0
    @test normal(p3) < 0
    @test centroid(p3) == Point(0.5, 0.5)
    @test Point(0.5, 0.5) ∈ p3
    @test Point(0.5, -0.5) ∉ p3
    @test Point(eps(), eps()) ∈ p3
    @test Point(-eps(), -eps()) ∉ p3

    p4 = ConvexPolygon(Point{3,Float64}.([(0,0,0), (1,0,0), (1,1,0), (0,1,0)]))
    @test area(p4) == 1.0
    @test normal(p4) == Point(0.0, 0.0, 1.0)
    @test centroid(p4) == Point(0.5, 0.5, 0.0)
    v4,conn4 = poly2mesh(p4)
    @test coordinates(p4)==v4
    @test conn4 == [1 2 3; 1 3 4]
    
    
    p5 = ConvexPolygon(Point{3,Float64}.( reverse([(0,0,0), (1,0,0), (1,1,0), (0,1,0)]) ))
    
    @test area(p5) == 1.0
    @test normal(p5) == Point(0.0, 0.0, -1.0)
    @test centroid(p4) == Point(0.5, 0.5, 0.0)
    
    
    
end
