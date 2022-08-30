export discrsurface
# Surface discretization

function discrsurface(tri, pts::AbstractVector{Point{3,Float64}})

    ntri = length(tri)

    npts = length(pts)

    # First of all get the bounding box of the triangles:
    xtmin = minimum(min(t[1][1], t[2][1], t[3][1]) for t in tri)
    xtmax = maximum(max(t[1][1], t[2][1], t[3][1]) for t in tri)

    ytmin = minimum(min(t[1][2], t[2][2], t[3][2]) for t in tri)
    ytmax = maximum(max(t[1][2], t[2][2], t[3][2]) for t in tri)

    ztmin = minimum(min(t[1][3], t[2][3], t[3][3]) for t in tri)
    ztmax = maximum(max(t[1][3], t[2][3], t[3][3]) for t in tri)

    # Same for the points
    xpmin = minimum(p[1] for p in pts)
    xpmax = maximum(p[1] for p in pts)

    ypmin = minimum(p[2] for p in pts)
    ypmax = maximum(p[2] for p in pts)

    zpmin = minimum(p[3] for p in pts)
    zpmax = maximum(p[3] for p in pts)


    xmin = min(xtmin, xpmin); ymin = min(ytmin, ypmin);
    zmin = min(ztmin, zpmin);
    
    xmax = max(xtmax, xpmax); ymax = min(ytmax, ypmax);
    zmax = min(ztmax, zpmax);
    
    Δ = 10*max(xmax-xmin, ymax-ymin, zmax-zmin)
    bbox = (x = (xmin-Δ, xmax+Δ), y = (ymin-Δ, ymax+Δ), z = (zmin-Δ, zmax+Δ))


    vor = voronoi3d(pts, bbox=bbox)
    # Chop each triangle with every polyhedron.

    # Each node of `pts` corresponds to a volume (polyhedron). We will
    # Get every triangle on the surface mesh and 
    trivor = Vector{TriangleFace{Point3{Float64}}}[]
    tidx = Vector{Int}[]
    for vol in vor.cells
        id = Int[]
        trim = TriangleFace{Point{3,Float64}}[]
        for i in 1:ntri
            t = tri[i]
            m = chopwithpolyhedron(vol, t)
            nm = length(m)
            if nm > 0
                for ti in m
                    if area(ti) > 0
                        push!(trim, ti)
                        push!(id, i)
                    end
                end
            end
        end
        push!(trivor, trim)
        push!(tidx, id)
    end

    return trivor, tidx
end


                      
