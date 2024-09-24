# Testing geometry transformations 

let
    # 2d
    p1 = SVec(0.0, 0.0)
    p2 = SVec(1.0, 0.0)
    p3 = SVec(1.0, 1.0)
    p4 = SVec(0.0, 1.0)
    
    @test translate(p1, SVec(0.,0.)) == p1
    @test translate(p1, SVec(1.,0.)) == p2
    @test translate(p1, SVec(1.,1.)) == p3
    
    # Testing reescaling points

    @test reescale(p3, scale=1/1000, origin=SVec(0.,0.), unit="mm", ounit="m") == p3
    @test reescale(p3, scale=1.0, origin=SVec(1.,1.), unit="mm", ounit="mm") == p1

    @test reescale(p3, 1.0, SVec(0., 0.)) == p3
    @test reescale(p3, 1.0, SVec(1., 1.)) == p1

    tri = Tri(p1, p2, p3)
    tri2 = reescale(tri, 1.0, SVec(0.0,0.0))
    
    # No transformation
    @test tri == tri2

    # Using scale and change of units
    tri2 = reescale(tri, scale=1/1000, origin=SVec(0.,0.),
                        unit="mm", ounit="m")
    @test tri2 == tri

    # Just changing the units
    tri2 = reescale(tri, scale=1.0, origin=SVec(0.,0.),
                        unit="mm", ounit="m")

    @test vertices(tri2) .* 1000 == vertices(tri)

    

    # Rotation matrix
    @test BuildingGeometry.rotationmatrix(π) == [-1.0  0.0;
                                                 0.0 -1.0]
    
    @test BuildingGeometry.rotationmatrix(π/2) ≈ [0.0  -1.0;
                                                  1.0 0.0]
    R = BuildingGeometry.rotationmatrix(π/4)
    @test rotate(p2*sqrt(2), R, SVec(0.,0.)) ≈ p3

    R = BuildingGeometry.rotationmatrix(1π)
    @test vertices(rotate(tri, R, SVec(0.,0.))) ≈ [p1, -p2, -p3]
    @test vertices(rotate(tri, 1π, SVec(0., 0.))) ≈ [p1, -p2, -p3]

    R = BuildingGeometry.rotationmatrix(π/2)
    @test vertices(rotate(tri, R, SVec(0.,0.))) ≈ [p1, rotate(p2,π/2,SVec(0.0,0.)),
                                                       rotate(p3,π/2,SVec(0.0,0.))]
    
    @test vertices(rotate(tri,π/2, SVec(0.,0.)))≈ [p1, rotate(p2,π/2,SVec(0.0,0.)),
                                                   rotate(p3,π/2,SVec(0.0,0.))]
    


end

let

    # 3d
    
    p1 = SVec(0.0, 0.0, 0.0)
    p2 = SVec(1.0, 0.0, 0.0)
    p3 = SVec(1.0, 1.0, 0.0)

    @test translate(p1, SVec(0.,0.,0.)) == p1
    @test translate(p1, SVec(1.,0.,0.)) == p2
    @test translate(p1, SVec(1.,1.,0.)) == p3
    
    # Testing reescaling points

    @test reescale(p3, scale=1/1000, origin=SVec(0.,0.,0.), unit="mm", ounit="m") == p3
    @test reescale(p3, scale=1.0, origin=SVec(1.,1.,0.), unit="mm", ounit="mm") == p1

    @test reescale(p3, 1.0, SVec(0., 0.,0.)) == p3
    @test reescale(p3, 1.0, SVec(1., 1.,0.)) == p1

    tri = Tri(p1, p2, p3)
    tri2 = reescale(tri, 1.0, SVec(0.0,0.0,0.))
    
    # No transformation
    @test tri == tri2

    # Using scale and change of units
    tri2 = reescale(tri, scale=1/1000, origin=SVec(0.,0.,0.),
                        unit="mm", ounit="m")
    @test tri2 == tri

    # Just changing the units
    tri2 = reescale(tri, scale=1.0, origin=SVec(0.,0.,0.),
                        unit="mm", ounit="m")

    @test vertices(tri2) .* 1000 == vertices(tri)

    
    w = SVec(0., 0., 1.)
    # Rotation matrix
    @test BuildingGeometry.rotationmatrix(1π,w) ≈ [-1.0  0.0  0.0
                                                   0.0 -1.0  0.0
                                                   0.0  0.0  1.0]
    
    @test BuildingGeometry.rotationmatrix(π/2,w) ≈ [0.0  -1.0 0.0
                                                  1.0 0.0  0.0
                                                  0.0 0.0  1.0]
    R = BuildingGeometry.rotationmatrix(π/4,w)
    @test rotate(p2*sqrt(2), R, SVec(0.,0.,0.)) ≈ p3

    R = BuildingGeometry.rotationmatrix(1π,w)
    @test vertices(rotate(tri, R, SVec(0.,0.,0.))) ≈ [p1, -p2, -p3]
    @test vertices(rotate(tri, 1π, w, SVec(0., 0.,0.))) ≈ [p1, -p2, -p3]
    zv = SVec(0.,0.,0.)
    R = BuildingGeometry.rotationmatrix(π/2,w)
    @test vertices(rotate(tri, R, SVec(0.,0.,0.))) ≈
        [p1, rotate(p2,π/2,SVec(0.0,0.,1.),zv),                                                       rotate(p3,π/2,SVec(0.0,0.,1.),zv)]
    
    @test vertices(rotate(tri,π/2, SVec(0.,0.,1.))) ≈
        [p1,rotate(p2,π/2,SVec(0.,0.,1.), zv),
         rotate(p3,π/2,SVec(0.,0.,1.), zv)]
    
    

    
end
