

let
    # Basic stuff. A single external surface
    pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(-1,1,0))]
    m =  buildsurface(tri, [(points=pts, tri=1:2, id=1:4)])
    Ne = length(pts)
    Ae = zeros(Ne)
    for n in m.nodes
        Ae[n.iex] += norm(n.A)
    end
    @test sum(Ae) ≈ 4.0
    @test all(Ae .≈ 1.0)


    ipts = Point.([(-0.75, 0,0), (0.75, 0,0)])
    m =  buildsurface(tri, [(points=pts, tri=1:2, id=1:4),
                            (points=ipts, tri=1:2, id=5:6)])
    Ne = length(pts)
    Ni = length(ipts)

    Ae = zeros(Ne)
    Ai = zeros(Ni)
    
    for n in m.nodes
        A = norm(n.A)
        Ae[n.iex] += A
        Ai[n.iin-4] += A
    end
    @test sum(Ae) ≈ 4.0
    @test sum(Ai) ≈ 4.0
    @test all(Ae .≈ 1.0)
    @test all(Ai .≈ 2.0)

    # Lets try something more complex
    tri = [Triangle((0,0,0), (1,0,0), (1,1,0)), Triangle((0,0,0), (1,1,0), (0,1,0))]
    Ne = 10
    Ni = 5
    epts = Point.(rand(Ne), rand(Ne), 0)
    ipts = Point.(rand(Ni), rand(Ni), 0)
    m =  buildsurface(tri, [(points=epts, tri=1:2, id=1:Ne),
                            (points=ipts, tri=1:2, id=1:Ni)])
    
    Ae = zeros(Ne)
    Ai = zeros(Ni)
    
    for n in m.nodes
        A = norm(n.A)
        Ae[n.iex] += A
        Ai[n.iin] += A
    end
    @test sum(Ae) ≈ 1.0
    @test sum(Ai) ≈ 1.0
    
    
    tri = [Triangle((0.,0.,0.), (1.,0.,0.), (1.,0.75,0.)),
              Triangle((0.,0.,0.), (1.,0.75,0.), (0.,0.75,0.)),
              Triangle((0.,0.75,0.), (1.,0.75,0.), (1.,1.,0.)),
              Triangle((0.,0.75,0.), (1.,1.,0.), (0.,1.,0.))]

    Ne = 50
    Ni = 20
    epts = Point.(rand(Ne), rand(Ne), 0)
    ipts = Point.(0.75*rand(Ni), 0.75*rand(Ni), 0)
    m =  buildsurface(tri, [(points=epts, tri=1:4, id=1:Ne),
                            (points=ipts, tri=1:2, id=1:Ni),
                            (tri=3:4, id=-1)])
    
    
    Ae = zeros(Ne)
    Ai = zeros(Ni)
    
    for n in m.nodes
        A = norm(n.A)
        Ae[n.iex] += A
        n.iin > 0 && (Ai[n.iin] += A)
    end
    @test sum(Ae) ≈ 1.0
    @test sum(Ai) ≈ 0.75
    
    pts = Point.([(-0.5,-0.5,0),(0.5,-0.5,0),(-0.5,0.5,0),(0.5,0.5,0)])
    tri = [Triangle((-1,-1,0), (1,-1,0), (1,1,0)), Triangle((-1,-1,0),(1,1,0),(-1,1,0))]

end

