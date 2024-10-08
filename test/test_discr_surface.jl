

function genpoints2d(W, H, nw, nh)
    Δh = H / nh
    Δw = W / nw
    
    h₁ = range(Δh/2, step=Δh, length=nh)
    w₁ = range(Δw/2, step=Δw, length=nw)
    return w₁, h₁
end
function genpoints(x,y, p0, u, v) 
    p = [(x1*u + y1*v) + p0 for x1 in x, y1 in y]
    return reshape(p, (length(x)*length(y),))
end


function gentripoints2d(W, H, nw, nh)
    
    x = range(zero(W), W, length=nw+1)
    y = range(zero(H), H, length=nh+1)
    
    return x,y 
    
end

function gentri(x,y, p0::SVec{Dim,T},
                u::SVec{Dim,T}, v::SVec{Dim,T}) where {Dim,T}
    
    nx = length(x)
    ny = length(y)
    
    np = nx*ny
    
    p = [(x1*u + y1*v + p0) for x1 in x, y1 in y]
    pidx = reshape(1:np, (nx, ny)) 
    
    pp = reshape(p, (nx*ny,))
    count = 1
    ntri = (nx-1)*(ny-1) * 2
    tri = Vector{Tri{Dim,T}}(undef, ntri)
    conn = zeros(Int, ntri, 3)
    for iy in 1:ny-1
        for ix in 1:nx-1
            tri[count] = Tri(p[ix,iy], p[ix+1,iy], p[ix+1,iy+1])
            tri[count+1] = Tri(p[ix,iy], p[ix+1,iy+1], p[ix,iy+1])
            conn[count,:] .= [pidx[ix,iy], pidx[ix+1,iy], pidx[ix+1,iy+1]]
            conn[count+1,:] .= [pidx[ix,iy], pidx[ix+1,iy+1], pidx[ix,iy+1]]
            count += 2
        end
    end
    
    return tri, pp, conn
end


connect_triangles(conn) = [connect( (r[1], r[2], r[3]) ) for r in eachrow(conn)]

let
    W  = 1.0
    H = 1.0
    At = W*H
    
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 1.0, 0.0)
    p0 = SVec(0.0, 0.0, 0.0)
    
    pxy = genpoints2d(W, H, 2, 2)
    pts = genpoints(pxy[1], pxy[2], p0, u, v)

    Ai = At / length(pts)
    
    txy = gentripoints2d(W, H, 1, 1)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    xtri, xidx = discrsurface(tri, pts)
    
    A = [sum(area.(t)) for t in xtri]
    
    # Test total area:
    @test sum(A) ≈ At
    @test all(A .≈ Ai)


    W  = 1.0
    H = 1.0
    At = W*H
    
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 1.0, 0.0)
    p0 = SVec(0.0, 0.0, 0.0)
    
    pxy = genpoints2d(W, H, 2, 2)
    pts = genpoints(pxy[1], pxy[2], p0, u, v)

    Ai = At / length(pts)
    
    txy = gentripoints2d(W, H, 10, 10)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    xtri, xidx = discrsurface(tri, pts)
    A = [sum(area.(t)) for t in xtri]
    @test sum(A) ≈ At
    @test all(A .≈ Ai)





    W  = 10.0
    H = 40.0
    At = W*H
    
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 0.0, 1.0)
    p0 = SVec(0.0, 0.0, 0.0)
    
    pxy = genpoints2d(W, H, 3, 4)
    pts = genpoints(pxy[1], pxy[2], p0, u, v)

    Ai = At / length(pts)
    
    txy = gentripoints2d(W, H, 2, 3)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    xtri, xidx = discrsurface(tri, pts)
    A = [sum(area.(t)) for t in xtri]
    @test sum(A) ≈ At
    @test all(A .≈ Ai)


    W  = 10.0
    H = 40.0
    At = W*H
    
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 0.0, 1.0)
    p0 = SVec(0.0, 0.0, 0.0)
    
    pxy = genpoints2d(W, H, 7, 13)
    pts = genpoints(pxy[1], pxy[2], p0, u, v)

    Ai = At / length(pts)
    
    txy = gentripoints2d(W, H, 13, 19)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    xtri, xidx = discrsurface(tri, pts)
    A = [sum(area.(t)) for t in xtri]
    @test sum(A) ≈ At
    @test all(A .≈ Ai)
    

    
    # Let's test wit a slightly more chanlengin case:
    # A square building
    
    W  = 10.0
    H = 40.0
    At = W*H*4
    
    pxy = genpoints2d(W, H, 7, 13)
    txy = gentripoints2d(W, H, 13, 19)

    u1 = SVec(1.0, 0.0, 0.0)
    v1 = SVec(0.0, 0.0, 1.0)
    p1 = SVec(0.0, 0.0, 0.0)
    pts1 = genpoints(pxy[1], pxy[2], p1, u1, v1)
    tri1, ptri1, conn1 = gentri(txy[1], txy[2], p1, u1, v1)

    u2 = SVec(0.0, 1.0, 0.0)
    v2 = SVec(0.0, 0.0, 1.0)
    p2 = SVec(W, 0.0, 0.0)
    pts2 = genpoints(pxy[1], pxy[2], p2, u2, v2)
    tri2, ptri2, conn2 = gentri(txy[1], txy[2], p2, u2, v2)
    
    u3 = SVec(-1.0, 0.0, 0.0)
    v3 = SVec(0.0, 0.0, 1.0)
    p3 = SVec(W, W, 0.0)
    pts3 = genpoints(pxy[1], pxy[2], p3, u3, v3)
    tri3, ptri3, conn3 = gentri(txy[1], txy[2], p3, u3, v3)
    
    
    u4 = SVec(0.0, -1.0, 0.0)
    v4 = SVec(0.0, 0.0, 1.0)
    p4 = SVec(0.0, W, 0.0)
    pts4 = genpoints(pxy[1], pxy[2], p4, u4, v4)
    tri4, ptri4, conn4 = gentri(txy[1], txy[2], p4, u4, v4)

    pts = [pts1; pts2; pts3; pts4]
    tri = [tri1; tri2; tri3; tri4]

    Ai = At / length(pts)

    xtri, xidx = discrsurface(tri, pts)
    A = [sum(area.(t)) for t in xtri]
    @test sum(A) ≈ At
    @test all(A .≈ Ai)

    # Now we will test the slicind stuff
    z = range(0.0, H, length=18)
    trilst, triidx = slicemesh(tri, z)

    A = [sum(area.(t)) for t in trilst]
    @test sum(A) ≈ At
    @test all(A .≈ z[2]*4*W)
    


    tri = [Tri(SVec(0.0,0,0), SVec(1.0,0,0), SVec(1.0,0,1)),
           Tri(SVec(0.0,0,0), SVec(1.0,0,1), SVec(0.0,0,1))]
    xt, xi = slicemesh(tri, [0.0, 0.5, 1.0])
    A = [sum(area.(t)) for t in xt]

    @test sum(A) ≈ 1.0
    @test all(A .≈ 0.5)

    z = range(0.0, 1.0, length=6)
    xt, xi = slicemesh(tri, z)
    A = [sum(area.(t)) for t in xt]
    @test sum(A) ≈ 1.0
    @test all(A .≈ 0.2)
    
    
    z = range(0.0, 1.0, length=3)
    txy = gentripoints2d(1.0, 1.0, 2, 4)
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 0.0, 1.0)
    p0 = SVec(0.0, 0.0, 0.0)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    trilst, triidx = slicemesh(tri, z)

    A = [sum(area.(t)) for t in trilst]
    @test sum(A) ≈ 1.0
    @test all(A .≈ z[2])


    z = range(0.0, 1.0, length=4)
    txy = gentripoints2d(1.0, 1.0, 2, 4)
    u = SVec(1.0, 0.0, 0.0)
    v = SVec(0.0, 0.0, 1.0)
    p0 = SVec(0.0, 0.0, 0.0)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    trilst, triidx = slicemesh(tri, z)

    A = [sum(area.(t)) for t in trilst]
    @test sum(A) ≈ 1.0
    @test all(A .≈ z[2])

end
