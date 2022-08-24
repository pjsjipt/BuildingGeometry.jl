import GeometryBasics: Point, Point2, Point3

export voronoi2d, voronoi3d,  VoronoiCell2d, makebbox


function voronoi2d(x,y; bbox=nothing, nd=8)

    length(x) == length(y) || error("Vectors x and y should have the same length")
    
    if isnothing(bbox)
        xmin,xmax = extrema(x)
        ymin,ymax = extrema(y)

        Δ = max(xmax-xmin, ymax-ymin)  # A length for the domain

        # Let's make a big bounding box
        bbox = (x=(xmin-nd*Δ,xmax+nd*Δ), y=(ymin-nd*Δ, ymax+nd*Δ))
    end

    #=
    The coordinates of the bounding box will make up 4 extra points
    at the beggining of the point list. This is done so that the
    volumes that contain the infinity point can be excluded.
    =#

    bx = repeat([bbox.x[1], bbox.x[2]], 2)
    by = repeat([bbox.y[1], bbox.y[2]], inner=2)
    nbb = length(bx) # Number of bounding box points
    px = [bx;x]
    py = [by;y]
    
    vor = pyvoronoi[]([px py])

    regions = Vector{Int}[]
    for idx in vor.point_region
        push!(regions, vor.regions[idx+1] .+ 1) # Python is 0-based indexing!
    end

    # Let's get the vertices
    vx = vor.vertices[:,1]
    vy = vor.vertices[:,2]

    # Let's get the ridges right
    ridges = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    # We only care about ridges that don't contain infinity points!
    for (p,v) in vor.ridge_dict
        if v[1] < 0 || v[2] < 0
            continue
        end
        k = (p[1]+1,p[2]+1) # Python indices...
        v = (v[1]+1,v[2]+1)
        ridges[k] = p
    end
    
    
    return vor, Point2.(px, py), Point2.(vx, vy), regions, ridges
    
end


function build_voronoi_mesh2d(nb, vor, pts, verts, regions, ridges)
    np = length(pts) - nb

    

end


struct VoronoiCell2d{Dim,T}
    p::Int
    verts::Vector{Int}
    edges::Vector{Tuple{Int,Int}}
    neighs::Vector{Int}
end

function bboxaux!(xlims, p)
    x = repeat([xlims[1], xlims[2]], inner=p[1], outer=p[2])
    p[1] = p[1] * 2
    p[2] = p[2] ÷ 2
    return x
end


function makebbox(ranges...) where {T}
    nr = length(ranges)
    p = [1,2^(nr-1)] # Outer and inner initial
    return map(x->bboxaux!(x, p), ranges)
    
end

function voronoi3d(x,y,z; bbox=nothing, nd=8)
    length(x) == length(y) == length(z) || error("Vectors x, y and z should have the same length")
    
    if isnothing(bbox)
        bbox = let
            xmin,xmax = extrema(x)
            ymin,ymax = extrema(y)
            zmin,zmax = extrema(z)
            Δ = max(xmax-xmin, ymax-ymin, zmax-zmin)  # A length for the domain
            # Let's make a big bounding box
            (x=(xmin-nd*Δ,xmax+nd*Δ), y=(ymin-nd*Δ, ymax+nd*Δ), z=(zmin-nd*Δ, zmax+nd*Δ))
        end
        
    end
    
    #=
    The coordinates of the bounding box will make up 4 extra points
    at the beggining of the point list. This is done so that the
    volumes that contain the infinity point can be excluded.
    =#
    #bx, by, bz = makebbox(bbox.x, bbox.y, bbox.z)
    xm = (bbox.x[1] + bbox.x[2])/2
    ym = (bbox.y[1] + bbox.y[2])/2
    zm = (bbox.z[1] + bbox.z[2])/2
    bx = [bbox.x[1], bbox.x[2], xm, xm, xm, xm]
    by = [ym, ym, bbox.y[1], bbox.y[2], ym, ym]
    bz = [zm, zm, zm, zm, bbox.z[1], bbox.z[2]]
    
    nbb = length(bx) # Number of bounding box points
    px = [bx;x]
    py = [by;y]
    pz = [bz;z]
    
    vor = pyvoronoi[]([px py pz])

    p = Point.(x,y,z)
    b = Point.(bx, by,bz)
    vx = vor.vertices[:,1]
    vy = vor.vertices[:,2]
    vz = vor.vertices[:,3]
    
    v = Point.(vx, vy, vz)

    pp = [b; p]
    regions = Vector{Int}[]
    for idx in vor.point_region
        push!(regions, vor.regions[idx+1] .+ 1) # Python is 0-based indexing!
    end

    # Let's get the ridges right
    ridges = Dict{Tuple{Int,Int}, Vector{Int}}()
    # We only care about ridges that don't contain infinity points!
    for (p,v) in vor.ridge_dict
        if minimum(v) < 0
            continue
        end
        k = (p[1]+1,p[2]+1) # Python indices...
        ridges[k] = v .+ 1
    end

    # Let's create the volumes

    # The ridges specify polygons between points
    plst = [Vector{Int}[] for i in 1:length(p)]
   
    for (pi, vi) in ridges
        i,k = pi
        if i > nbb
            push!(plst[i-nbb], vi)
        end
        if k > nbb
            push!(plst[k-nbb], vi)
        end
    end

    #faces = [[ConvexPolygon(v[idx]) for idx in ff] for ff in plst]

    conn = volume2mesh.(plst)

    return vor, p,b,v,pp,regions, ridges, plst, conn
    
    
    
end

function volume2mesh(faceidx)

    # Convert to triangles
    nf = length(faceidx)
    ntri = sum(length.(faceidx)) - 2*nf

    conn = zeros(Int, ntri, 3)

    idx = 1
    for face in faceidx
        # Each face is a convex polygon
        nt = length(face)
        for k in 2:nt-1
            conn[idx,1] = face[1]
            conn[idx,2] = face[k]
            conn[idx,3] = face[k+1]
            idx += 1
        end
    end

    return conn
end
