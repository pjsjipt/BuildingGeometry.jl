export discrsurface
# Surface discretization

function discrsurface(tri, pts::AbstractVector{Point{3,Float64}};
                      rtol=1e-8, bbox=nothing, nd=8)

    ntri = length(tri)
    npts = length(pts)

    # First of all get the bounding box of the triangles:
    bbpts = boundingbox(pts)
    if isnothing(bbox)
        bbtri = boundingbox(tri)
        bbox1 = boundingbox([bbtri.min, bbtri.max, bbpts.min, bbpts.max])
        Δ = norm(bbox1.max-bbox1.min)
        u = Vec(Δ, Δ, Δ)
        bbox = Box(bbox1.min - nd*u, bbox1.max + nd*u)
    end

    Lref = norm(bbpts.max-bbpts.min)
    
    atol = rtol * Lref/2

    vor = voronoi3d(pts, bbox=bbox)
    # Chop each triangle with every polyhedron.

    # Each node of `pts` corresponds to a volume (polyhedron). We will
    # Get every triangle on the surface mesh and 
    trivor = Vector{Triangle{3,Float64}}[]
    tidx = Vector{Int}[]
    for vol in vor.cells
        id = Int[]
        trim = Triangle{3,Float64}[]
        for (i,t) in enumerate(tri)
            t = tri[i]
            m = chopwithpolyhedron(vol, t, atol=atol)
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


                      
