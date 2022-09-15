

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
    addforcecontrib!(Fe, m, (1,2,3,4,5,6); sgn=1, side=1)

    @test all(Fe[1,:] .≈ 0)  # No forces in x
    @test all(Fe[2,:] .≈ 0)  # No forces in y
    @test all(Fe[3,1:Ne] .≈ -1/Ne)

    # Let's test the moments
    @test all(Fe[6,:] .≈ 0)  # No moments in z
    Ae = 1/Ne
    @test all(Fe[4,1:Ne] .≈ -Ae .* y)
    @test all(Fe[5,1:Ne] .≈ +Ae .* x)

 
    Fi = zeros(6,10)
    addforcecontrib!(Fi, m, (1,2,3,4,5,6); sgn=1, side=2)
    @test all(Fi[1,:] .≈ 0)  # No forces in x
    @test all(Fi[2,:] .≈ 0)  # No forces in y
    @test all(Fi[3,Ne .+ (1:Ni)] .≈ -1/Ni)

    # Let's test the moments
    @test all(Fi[6,:] .≈ 0)  # No moments in z
    Ai = 1/Ni
    @test all(Fi[4,Ne.+(1:Ni)] .≈ -Ai .* [0.5, 0.5])
    @test all(Fi[5,Ne.+(1:Ni)] .≈ +Ai .* [0.25, 0.75])
    
    # Let's repeat but now the surface is on plane x-z
    
    epts = Point.(x, 0.0, y); Ne = length(epts)
    ipts = Point.([(0.2,0.0,0.0), (0.8,0.0,0.0)]); Ni = length(ipts)
    
    tri = [Triangle((0,0,0), (1,0,0), (1,0,1)),
           Triangle((0,0,0), (1,0,1), (0,0,1))]

    m = buildsurface(tri, [(points=epts,tri=1:2, id=1:Ne),
                           (points=ipts,tri=1:2, id=Ne .+ (1:Ni))])

    Fe = zeros(6,10)
    addforcecontrib!(Fe, m, (1,2,3,4,5,6); sgn=1, side=1)

    @test all(Fe[1,1:Ne] .≈ 0)  # No forces in x
    @test all(Fe[2,1:Ne] .≈ +1/Ne)  # No forces in y
    @test all(Fe[3,1:Ne] .≈ 0)
    Ae = 1/Ne
    Ai = 1/Ni
    @test all(Fe[4,1:Ne] .≈ -y * Ae) 
    @test all(Fe[5,1:Ne] .≈ 0)
    @test all(Fe[6,1:Ne] .≈ +x * Ae) 
   
    Fi = zeros(6,10)
    addforcecontrib!(Fi, m, (1,2,3,4,5,6); sgn=1, side=2)
     
    @test all(Fi[1,Ne .+ (1:Ni)] .≈ 0)  # No forces in x
    @test all(Fi[2,Ne .+ (1:Ni)] .≈ [+0.5, +0.5])  # 
    @test all(Fi[3,Ne .+ (1:Ni)] .≈ 0)  # No forces in z
    @test all(Fi[4,Ne .+ (1:Ni)] .≈ -0.5*0.5)  # No forces in x
    @test all(Fi[5,Ne .+ (1:Ni)] .≈ 0)  # No forces in x
    @test all(Fi[6,Ne .+ (1:Ni)] .≈ [+0.125, +0.375])  # No forces in x

    # Now I will try to generate a building with floor floors.
    # The building will have 4 sides, 8m wide and 16m high
    # Each face has 8 external nodes and 4 internal nodes
    # External nodes
    epts = Point.([(-2,-4,2),(2,-4,2), #F1
                   (-2,-4,6),(2,-4,6), #F1
                   (-2,-4,10),(2,-4,10), #F1
                   (-2,-4,14),(2,-4,14), #F1
                   (4,-2,2),(4,2,2), #F2
                   (4,-2,6),(4,2,6), #F2
                   (4,-2,10),(4,2,10), #F2
                   (4,-2,14),(4,2,14), #F2
                   (-2,4,2),(2,4,2), #F3
                   (-2,4,6),(2,4,6), #F3
                   (-2,4,10),(2,4,10), #F3
                   (-2,4,14),(2,4,14), #F3
                   (-4,-2,2),(-4,2,2), #F4
                   (-4,-2,6),(-4,2,6), #F4
                   (-4,-2,10),(-4,2,10), #F4
                   (-4,-2,14),(-4,2,14)]) #F4
    # Internal nodes
    ipts = Point.([(-2,-4,4),(2,-4,4), #F1
                   (-2,-4,12),(2,-4,12), #F1
                   (4,-2,4),(4,2,4), #F2
                   (4,-2,12),(4,2,12), #F2
                   (-2,4,4),(2,4,4), #F3
                   (-2,4,12),(2,4,12), #F3
                   (-4,-2,4),(-4,2,4),#F4
                   (-4,-2,12),(-4,2,12)])#F4

    # Indexes of external nodes on each face
    eidx = [1:8, 9:16, 17:24, 25:32]
    # Indexes of internal nodes on each face
    iidx = [1:4, 5:8, 9:12, 13:16]

    # Let's get the faces:
    tri = [Triangle((-4,-4,0),(4,-4,0),(4,-4,16)), # F1
           Triangle((-4,-4,0),(4,-4,16),(-4,-4,16)), # F1
           Triangle((4,-4,0),(4,4,0),(4,4,16)), # F2
           Triangle((4,-4,0),(4,4,16),(4,-4,16)), # F2
           Triangle((4,4,0),(-4,4,0),(-4,4,16)), # F3
           Triangle((4,4,0),(-4,4,16),(4,4,16)), # F3
           Triangle((-4,4,0),(-4,-4,0),(-4,-4,16)), # F4
           Triangle((-4,4,0),(-4,-4,16),(-4,4,16))] # F4
    triidx = [1:2, 3:4, 5:6, 7:8]
    
    Ne = 8*4
    Ni = 4*4
    Nt = Ne + Ni
                             
    faces = [buildsurface(tri, [(points=epts[eidx[i]], tri=triidx[i], id=eidx[i]),
                                (points=ipts[iidx[i]], tri=triidx[i], id=iidx[i].+Ne)])
             for i in 1:4]

    bmsh = mergemeshes(faces)

    Ae = zeros(Ne)
    Ai = zeros(Ni)

    # Let's calculate the areas
    for n in bmsh.nodes
        e = nodeside(n,1)
        i = nodeside(n,2)
        An = norm(nodearea(n))
        Ae[e] += An
        Ai[i-Ne] += An
    end

    @test all(Ae .≈ 16)
    @test all(Ai .≈ 32)

    # I will now slice the buiding
    zh = [0.0, 4.0, 8.0, 12.0, 16.0]
    floors = buildingslice(bmsh, zh)

    # Let's check the areas of each node in each slice
    nz = length(zh)-1
    A = zeros(Nt,nz)
    for (k,f) in enumerate(floors)
        a = norm.(nodearea.(f.nodes))
        e = nodeside.(f.nodes, 1)
        i = nodeside.(f.nodes, 2)
        for j in eachindex(a)
            A[e[j],k] += a[j]
            A[i[j],k] += a[j]
        end
    end

    for i in 1:4
        @test sum(A[:,i]) ≈ 8*4 * 4 * 2
    end
    # Each cell should have an area = 16 or 0
    @test sum(A.≈16) == 32*2

    # Agora, vamos trabalhar com as forças
    MF1 = forcematrix(Nt, bmsh, (1,2,3,4,5,6); sgn=1, side=1) 
    MF2 = forcematrix(Nt, bmsh, (1,2,3,4,5,6); sgn=1, side=2)
    @test all(MF1[:,Ne .+ (1:Ni)] .== 0)  # Internal nodes should not contribute
    @test all(MF2[:,1:Ne] .== 0) # External nodes should not contribute
    # For constant pressure, all forces should add up to zero
    @test all(isapprox.(sum(MF1,dims=2), 0, atol=1e-8))
    @test all(isapprox.(sum(MF2,dims=2), 0, atol=1e-8))

    # Let's check per face
    # Face 1. Fx = 0, Fy = 128, Fz = 0, Mx = 128*8, My=Mz=0
    idx = eidx[1]
    @test isapprox(sum(MF1[1,idx]),0, atol=1e-8)
    @test sum(MF1[2,idx]) ≈ +8*16
    @test isapprox(sum(MF1[3,idx]),0, atol=1e-8)
    @test sum(MF1[4,idx]) ≈ -8*16 * 8
    @test isapprox(sum(MF1[5,idx]),0, atol=1e-8)
    @test isapprox(sum(MF1[6,idx]),0, atol=1e-8)

    # Face 2. Fx = 128, Fy = 0
    idx = eidx[2]
    @test sum(MF1[1,idx]) ≈ -8*16
    @test isapprox(sum(MF1[2,idx]),0, atol=1e-8)
    @test isapprox(sum(MF1[3,idx]),0, atol=1e-8)
    @test isapprox(sum(MF1[4,idx]),0, atol=1e-8)
    @test sum(MF1[5,idx]) ≈ -8*16 * 8
    @test isapprox(sum(MF1[6,idx]),0, atol=1e-8)

    # Let's test  FM2, face 1
    idx = iidx[1] .+ Ne
    @test isapprox(sum(MF2[1,idx]),0, atol=1e-8)
    @test sum(MF2[2,idx]) ≈ +8*16
    @test isapprox(sum(MF2[3,idx]),0, atol=1e-8)
    @test sum(MF2[4,idx]) ≈ -8*16 * 8
    @test isapprox(sum(MF2[5,idx]),0, atol=1e-8)
    @test isapprox(sum(MF2[6,idx]),0, atol=1e-8)


    # Let's try to do the samething for slices.
    # We will calculate for each slice and rebuild the MF1 and MF2 matrices
    SF1 = forcematrix(Nt, floors, (1,2,3,4,5,6), side=1)
    SF2 = forcematrix(Nt, floors, (1,2,3,4,5,6), side=2)
    kf = [(1:4) .+ i for i in 0:4:23 ]
    SF1i = vcat([sum(SF1[k,:],dims=1) for k in kf]...)
    SF2i = vcat([sum(SF2[k,:],dims=1) for k in kf]...)
    @test maximum(abs, SF1i-MF1) < 1e-8
    @test maximum(abs, SF2i-MF2) < 1e-8
    # Let-'s test Fy and Mx for second floor first face, epts[3]
    @test SF1[6,3] ≈ +16
    @test SF1[14,3] ≈ -16*6
    
    
    SF1 = forcematrix(Nt, floors, (1,2,3,4,5,6), side=1, interleaved=true)
    SF2 = forcematrix(Nt, floors, (1,2,3,4,5,6), side=2, interleaved=true)
    kf = [(1:6:24) .+ i for i in 0:5 ]
    SF1i = vcat([sum(SF1[k,:],dims=1) for k in kf]...)
    SF2i = vcat([sum(SF2[k,:],dims=1) for k in kf]...)
    @test maximum(abs, SF1i-MF1) < 1e-8
    @test maximum(abs, SF2i-MF2) < 1e-8
    # Let-'s test the Fy for second floor, epts[3]
    @test SF1[8,3] ≈ +16
    @test SF1[10,3] ≈ -16*6
    
    
    
    
end
