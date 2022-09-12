

# Testing forces stuff


let
    x1 = [0.25, 0.75]
    y1 = [0.125, 0.375, 0.625, 0.875]

    x = repeat(x1, 4)
    y = repeat(y1, inner=2)
    
    epts = Point.(x, y, 0.0); Ne = length(epts)
    ipts = Point.([(0.2,0.5,0.0), (0.8,0.5,0.0)]); Ni = length(ipts)
    
    tri = [Triangle((0,0,0), (1,0,0), (1,1,0)),
           Triangle((0,0,0), (1,1,0), (0,1,0))]

    m = buildsurface(tri, [(points=epts,tri=1:2, id=1:Ne),
                           (points=ipts,tri=1:2, id=Ne .+ (1:Ni))])

    Ae = zeros(Ne)
    Ai = zeros(Ni)

    for node in m.nodes
        A = norm(nodearea(node))
        e = nodeside(node,1)
        i = nodeside(node,2)
        Ae[e] += A
        Ai[i-Ne] += A
    end

    @test sum(Ae) ≈ 1.0
    @test sum(Ai) ≈ 1.0

    @test all(Ae .≈ 1/Ne)
    @test all(Ai .≈ 1/Ni)

    Fe = zeros(6,10)
    addforcecontrib!(Fe, m.nodes, (1,2,3,4,5,6); sgn=1, side=1)

    @test all(Fe[1,:] .≈ 0)  # No forces in x
    @test all(Fe[2,:] .≈ 0)  # No forces in y
    @test all(Fe[3,1:Ne] .≈ 1/Ne)

    # Let's test the moments
    @test all(Fe[6,:] .≈ 0)  # No moments in z
    Ae = 1/Ne
    @test all(Fe[4,1:Ne] .≈ Ae .* y)
    @test all(Fe[5,1:Ne] .≈ -Ae .* x)

 
    Fi = zeros(6,10)
    addforcecontrib!(Fi, m.nodes, (1,2,3,4,5,6); sgn=1, side=2)
    @test all(Fi[1,:] .≈ 0)  # No forces in x
    @test all(Fi[2,:] .≈ 0)  # No forces in y
    @test all(Fi[3,Ne .+ (1:Ni)] .≈ 1/Ni)

    # Let's test the moments
    @test all(Fi[6,:] .≈ 0)  # No moments in z
    Ai = 1/Ni
    @test all(Fi[4,Ne.+(1:Ni)] .≈ Ai .* [0.5, 0.5])
    @test all(Fi[5,Ne.+(1:Ni)] .≈ -Ai .* [0.25, 0.75])
    
    # Let's repeat but now the surface is on plane x-z
    
    epts = Point.(x, 0.0, y); Ne = length(epts)
    ipts = Point.([(0.2,0.0,0.0), (0.8,0.0,0.0)]); Ni = length(ipts)
    
    tri = [Triangle((0,0,0), (1,0,0), (1,0,1)),
           Triangle((0,0,0), (1,0,1), (0,0,1))]

    m = buildsurface(tri, [(points=epts,tri=1:2, id=1:Ne),
                           (points=ipts,tri=1:2, id=Ne .+ (1:Ni))])

    Fe = zeros(6,10)
    addforcecontrib!(Fe, m.nodes, (1,2,3,4,5,6); sgn=1, side=1)

    @test all(Fe[1,1:Ne] .≈ 0)  # No forces in x
    @test all(Fe[2,1:Ne] .≈ -1/Ne)  # No forces in y
    @test all(Fe[3,1:Ne] .≈ 0)
    Ae = 1/Ne
    Ai = 1/Ni
    @test all(Fe[4,1:Ne] .≈ y * Ae) 
    @test all(Fe[5,1:Ne] .≈ 0)
    @test all(Fe[6,1:Ne] .≈ -x * Ae) 
   
    Fi = zeros(6,10)
    addforcecontrib!(Fi, m.nodes, (1,2,3,4,5,6); sgn=1, side=2)
     
    @test all(Fi[1,Ne .+ (1:Ni)] .≈ 0)  # No forces in x
    @test all(Fi[2,Ne .+ (1:Ni)] .≈ [-0.5, -0.5])  # 
    @test all(Fi[3,Ne .+ (1:Ni)] .≈ 0)  # No forces in z
    @test all(Fi[4,Ne .+ (1:Ni)] .≈ 0.5*0.5)  # No forces in x
    @test all(Fi[5,Ne .+ (1:Ni)] .≈ 0)  # No forces in x
    @test all(Fi[6,Ne .+ (1:Ni)] .≈ [-0.125, -0.375])  # No forces in x
   
end
