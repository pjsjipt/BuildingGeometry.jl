
let

    # Let's create a very simple geometry

    tri =[Tri(SVec(0.0,0.0,0.0),
              SVec(1.0,0.0,0.0), 
              SVec(1.0,1.0,0.0)),
          Tri(SVec(0.0,0.0,0.0),
              SVec(1.0,1.0,0.0),
              SVec(0.0,1.0,0.0)),
          Tri(SVec(0.0,1.0,0.0),
              SVec(1.0,1.0,0.0), 
              SVec(1.0,2.0,0.0)),
          Tri(SVec(0.0,0.0,0.0),
              SVec(1.0,2.0,0.0),
              SVec(0.0,2.0,0.0))]

    cpts = centroid.(tri)

    A = normalarea.(tri)

    msh = [nodeinfo(A[1], cpts[1], (1,-1)),
           nodeinfo(A[1], cpts[2], (2,-1)),
           nodeinfo(A[1], cpts[3], (3,5)),
           nodeinfo(A[1], cpts[4], (4,5))]
    
    sdata = PressToMesh(msh, 1)

    @test sdata.A == ones(4,1)
    @test sdata.idx == hcat(1:4)

    # Now we test the application
    cp1 = randn(10)

    @test all(cp1[1:4] .== assemble(sdata,cp1))
    
    # Let's check this for time series
    cp2 = randn(100,10)
    @test all(cp2[:,1:4] .== assemble(sdata,cp2))


    # Let's check the other side
    sidata = PressToMesh(msh,2)
    @test sidata.A[:,1] == [0.0, 0.0, 1.0, 1.0]
    @test sidata.idx[:,1] == [1,1,5,5]

    @test all([0.0, 0.0, cp1[5], cp1[5]] .== assemble(sidata, cp1))

    @test all([zeros(100) zeros(100) cp2[:,5] cp2[:,5]] .== assemble(sidata,cp2))
end
