

let

    #pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    #tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(0,1,0))]
    #vor = discrsurface(tri, 1:2, pts)

    # Triangle intersection 
    # First trivial case: no intersection whatsoever
    tri1 = Triangle((0,0,0),(1,0,0),(1,1,0))
    tri2 = Triangle((10,0,0),(11,0,0),(11,1,0))
    pts = intersect_tri(tri1, tri2)

    @test length(pts) == 0

    # Second trivial test: one triangle is inside the other
    tri1 = Triangle((0,0,0),(1,0,0),(1,1,0))
    tri2 = Triangle((-2,-1,0),(2,-1,0), (2,3,0))
    pts = intersect_tri(tri1, tri2)
    
    @test length(pts) == 3

    idx = findfirst(isequal(Point(0,0,0)), pts)
    @test !isnothing(idx)
    ifun(i,n=3) = (i-1) % n + 1

    v = vertices(tri1)
    for i in 1:3
        @test v[i] == pts[ifun(idx + i - 1, 3)]
    end
    
    # Maybe we should tilt this
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((-2,-1,-3),(2,-1,1), (2,3,5))
    pts = intersect_tri(tri1, tri2)
    
    @test length(pts) == 3

    idx = findfirst(isequal(Point(0,0,0)), pts)
    @test !isnothing(idx)

    v = vertices(tri1)
    for i in 1:3
        @test v[i] == pts[ifun(idx + i - 1, 3)]
    end

    # Let's invert the search
    pts = intersect_tri(tri2, tri1)
    
    @test length(pts) == 3

    idx = findfirst(isequal(Point(0,0,0)), pts)
    @test !isnothing(idx)

    v = vertices(tri1)
    for i in 1:3
        @test v[i] == pts[ifun(idx + i - 1, 3)]
    end


    # Lets see what happens if the triangles share part of an edge
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((0,0,0),(2,0,2),(2,2,4))

    pts = intersect_tri(tri1, tri2)

    @test length(pts) == 3

    idx = findfirst(isequal(Point(0,0,0)), pts)
    @test !isnothing(idx)
    v = vertices(tri1)
    for i in 1:3
        @test v[i] == pts[ifun(idx + i - 1, 3)]
    end

    # Let's try something harder
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((0.5,-0.5,0.),(1.5,0.,1.5),(0.5, 1.0, 1.5))

    pts = intersect_tri(tri1, tri2)

    @test length(pts) == 5
    idx = findfirst(isapprox(Point(1,0,1)), pts)
    @test pts[idx] ≈ Point(1,0,1)
    @test pts[ifun(idx+1,5)] ≈ Point(1.0, 0.5, 1.5)
    @test pts[ifun(idx+2,5)] ≈ Point(0.75, 0.75, 1.5)
    @test pts[ifun(idx+3,5)] ≈ Point(0.5, 0.5, 1.0)
    @test pts[ifun(idx+4,5)] ≈ Point(0.5, 0.0, 0.5)

    # We will invert this again...
    pts = intersect_tri(tri2, tri1)

    @test length(pts) == 5
    idx = findfirst(isapprox(Point(1,0,1)), pts)
    @test pts[idx] ≈ Point(1,0,1)
    @test pts[ifun(idx+1,5)] ≈ Point(1.0, 0.5, 1.5)
    @test pts[ifun(idx+2,5)] ≈ Point(0.75, 0.75, 1.5)
    @test pts[ifun(idx+3,5)] ≈ Point(0.5, 0.5, 1.0)
    @test pts[ifun(idx+4,5)] ≈ Point(0.5, 0.0, 0.5)


    #===================================================================
    # Testin mesh intersection
    ====================================================================#

    pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(0,1,0))]
    tvor, tidx = discrsurface(tri, 1:2, pts)

end

