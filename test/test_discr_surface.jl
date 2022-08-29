

let
    function genpoints2d(W, H, nw, nh)
        Δh = H / nh
        Δw = W / nw

        h₁ = range(Δh/2, step=Δh, length=nh)
        w₁ = range(Δw/2, step=Δw, length=nw)
        return w₁, h₁
    end
    function genpoints(x,y, p0, u, v) 
        p = [x1*u + y1*v + p0 for x1 in x, y1 in y]
        return reshape(p, (length(x)*length(y),))
    end
    
                       
    function gentripoints2d(W, H, nw, nh)

        Δx = W / nw
        Δy = H / nw

        x = range(0.0, W, length=nw+1)
        y = range(0.0, W, length=nh+1)

        return x,y 
        
    end

    function gentri(x,y, p0, u, v)

        nx = length(x)
        ny = length(y)

        np = nx*ny

        p = [(x1*u + y1*v + p0) for x1 in x, y1 in y]
        pidx = reshape(1:np, (nx, ny)) 
       
        pts = reshape(p, (nx*ny,))
        count = 1
        ntri = (nx-1)*(ny-1) * 2
        tri = Vector{TriangleFace{Point{3,Float64}}}(undef, ntri)
        conn = zeros(Int, ntri, 3)
        for iy in 1:ny-1
            for ix in 1:nx-1
                tri[count] = TriangleFace(p[ix,iy], p[ix+1,iy], p[ix+1,iy+1])
                tri[count+1] = TriangleFace(p[ix,iy], p[ix+1,iy+1], p[ix,iy+1])
                conn[count,:] .= [pidx[ix,iy], pidx[ix+1,iy], pidx[ix+1,iy+1]]
                conn[count+1,:] .= [pidx[ix,iy], pidx[ix+1,iy+1], pidx[ix,iy+1]]
                count += 2
            end
        end

        return tri, pts, conn
    end

    W  = 1.0
    H = 1.0

    u = Point3(1.0, 0.0, 0.0)
    v = Point3(0.0, 1.0, 0.0)
    p0 = Point3(0.0, 0.0, 0.0)
    
    pxy = genpoints2d(W, H, 2, 2)
    pts = genpoints(pxy[1], pxy[2], p0, u, v)

    txy = gentripoints2d(W, H, 1, 1)
    tri, ptri, conn = gentri(txy[1], txy[2], p0, u, v)
    
    #
        

    
        

end
