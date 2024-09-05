import LinearAlgebra: norm
import StaticArrays: SVector

cmppoint(p1, p2) = Meshes.isapproxzero(norm(p1-p2))
let

    #pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    #tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(0,1,0))]
    #vor = discrsurface(tri, 1:2, pts)

    # Triangle intersection 
    # First trivial case: no intersection whatsoever
    tri1 = Triangle((0,0,0),(1,0,0),(1,1,0))
    tri2 = Triangle((10,0,0),(11,0,0),(11,1,0))
    tris = intersect_tri(tri1, tri2)

    @test length(tris) == 0

    # Second trivial test: one triangle is inside the other
    tri1 = Triangle((0,0,0),(1,0,0),(1,1,0))
    tri2 = Triangle((-2,-1,0),(2,-1,0), (2,3,0))
    tris = intersect_tri(tri1, tri2)
    
    @test length(tris) == 1

    pts = vertices(tris[begin])
    idx = findfirst(p->cmppoint(p,Point(0,0,0)), pts)
    @test !isnothing(idx)
    ifun(i,n=3) = (i-1) % n + 1

    v = vertices(tri1)
    for i in 1:3
        @test cmppoint(v[i], pts[ifun(idx + i - 1, 3)])
    end
    
    # Maybe we should tilt this
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((-2,-1,-3),(2,-1,1), (2,3,5))
    tris = intersect_tri(tri1, tri2)
    
    @test length(tris) == 1
    pts = vertices(tris[begin])
    idx = findfirst(p -> cmppoint(p,Point(0,0,0)), pts)
    @test !isnothing(idx)

    v = vertices(tri1)
    for i in 1:3
        @test cmppoint(v[i], pts[ifun(idx + i - 1, 3)])
    end

    # Let's invert the search
    tris = intersect_tri(tri2, tri1)
    
    @test length(tris) == 1
    pts = vertices(tris[begin])
    idx = findfirst(p->cmppoint(p,Point(0,0,0)), pts)
    @test !isnothing(idx)

    v = vertices(tri1)
    for i in 1:3
        @test cmppoint(v[i], pts[ifun(idx + i - 1, 3)])
    end


    # Lets see what happens if the triangles share part of an edge
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((0,0,0),(2,0,2),(2,2,4))

    tris = intersect_tri(tri1, tri2)

    @test length(tris) == 1
    pts = vertices(tris[begin])
    idx = findfirst(p->cmppoint(p,Point(0,0,0)), pts)
    @test !isnothing(idx)
    v = vertices(tri1)
    for i in 1:3
        @test cmppoint(v[i], pts[ifun(idx + i - 1, 3)])
    end

    # Let's try something harder
    tri1 = Triangle((0,0,0),(1,0,1),(1,1,2))
    tri2 = Triangle((0.5,-0.5,0.),(1.5,0.,1.5),(0.5, 1.0, 1.5))

    tris = intersect_tri(tri1, tri2)

    px = [Point(1,0,1), Point(1.0, 0.5, 1.5), Point(0.75, 0.75, 1.5),
          Point(0.5, 0.5, 1.0), Point(0.5, 0.0, 0.5)]

    poly = ConvexPolygon(px)

    @test area(poly) ≈ sum(area.(tris))

    # We will invert this again...
    tris = intersect_tri(tri2, tri1)
    @test area(poly) ≈ sum(area.(tris))
    


    #===================================================================
    # Testing mesh intersection
    ====================================================================#

    pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(-1,1,0))]
    evor, eidx = discrsurface(tri, 1:2, pts)

    ptsi = Point.([(-0.6,0.0,0.0), (0.6,0.0,0.0)])
    ivor, iidx = discrsurface(tri, 1:2, ptsi)
    
    msh = eltype(tri)[]
    nodes = NodeInfo{Float64,Tuple{Int,Int}}[]

    Ne = length(pts)
    Ni = length(ptsi)

    for e in 1:Ne
        for i in 1:Ni
            intersectmesh!(msh, nodes, e, evor[e], eidx[e],
                           i, ivor[i], iidx[i])
        end
    end
    

    # Let's test the basics
    @test length(evor) == length(pts) == length(eidx)
    @test length(ivor) == length(ptsi) == length(iidx)
    
    # Let's check the areas:
    Ae = zeros(Ne)
    Ai = zeros(Ni)

    for node in nodes
        A = norm(node.A)
        Ae[node.side[1]] += A
        Ai[node.side[2]] += A
    end

    @test sum(Ae) ≈ 4.0
    @test sum(Ai) ≈ 4.0
    @test all(Ae .≈ 1.0)
    @test all(Ai .≈ 2.0)
    
    
    # Let's get more points and a more complex geometry

    tri = [Triangle((0,0,0), (1,0,0), (1,1,0)), Triangle((0,0,0),(1,1,0),(0,1,0))]
    epts = Point.(rand(10), rand(10), 0.0)
    ipts = Point.(rand(5), rand(5), 0.0)

    evor, eidx = discrsurface(tri, 1:2, epts)
    ivor, iidx = discrsurface(tri, 1:2, ipts)

    TriFace = eltype(tri)

    msh = eltype(tri)[]
    nodes = NodeInfo{Float64,Tuple{Int,Int}}[]

    Ne = length(epts)
    Ni = length(ipts)

    for e in 1:Ne
        for i in 1:Ni
            intersectmesh!(msh, nodes, e, evor[e], eidx[e],
                           i, ivor[i], iidx[i])
        end
    end
    
    # Let's check the areas:
    Ae = zeros(Ne)
    Ai = zeros(Ni)

    for node in nodes
        A = norm(node.A)
        Ae[node.side[1]] += A
        Ai[node.side[2]] += A
    end

    Ave = [sum(area.(ev)) for ev in evor] 
    Avi = [sum(area.(iv)) for iv in ivor]
    
    @test sum(Ave) ≈ 1.0m²
    @test sum(Avi) ≈ 1.0m²
    @test sum(Ae) ≈ 1.0m²
    @test sum(Ai) ≈ 1.0m²
    
    @test Ave ≈ Ae
    @test Avi ≈ Ai
    

end

