
import GeometryBasics: AbstractPoint, Point, Point3

export voronoi2d, voronoi3d,  VoronoiCell2d, makebbox
export getpoint, getvertex


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
    conn = [Int[] for i in 1:length(p)]
    
    for (pi, vi) in ridges
        i,k = pi
        i1 = i > nbb ? i-nbb : 0
        k1 = k > nbb ? k-nbb : 0
        
        # Now we will create a polygon and check its orientation
        # with respect to the points. If the normal of the polygon
        # points towards point `i` it should be reversed and for
        # point `k` the orientation is correct.
        pface = ConvexPolygon(v[vi])
        n⃗ = normal(pface)
        # Let's create a vector point pointing from `i` to any vertex
        # of the polygonal face:
        u⃗ = v[vi[1]] - pp[i]
        α = u⃗ ⋅ n⃗
        if i > nbb
            if α > 0
                push!(plst[i1], vi)
            else
                push!(plst[i1], reverse(vi))
            end
            push!(conn[i1], k1)
            
        end
        if k > nbb
            if α > 0
                push!(plst[k-nbb], reverse(vi))
            else
                push!(plst[k-nbb], vi)
            end
            push!(conn[k1], i1)
        end
    end

    # Setup data structure
    cells = [ConvexPolyhedron(v, iface) for iface in plst]

    return VoronoiMesh(p, v, cells, conn)
    
    
    
end

function voronoi3d(pts::AbstractVector{<:AbstractPoint{3}}; bbox=nothing, nd=8)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    z = [p[3] for p in pts]
    return voronoi3d(x, y, z; bbox=bbox, nd=nd)
end

             
    
    
struct VoronoiMesh{Dim,T,P<:AbstractPoint{Dim,T},
                   VP<:AbstractVector{P},
                   VV<:AbstractVector{P},CP,
                   VPO<:AbstractVector{CP}}
    points::VP
    vertices::VV
    cells::VPO
    conn::Vector{Vector{Int}}
end

function Base.show(io::IO, msh::VoronoiMesh{Dim,T}) where {Dim,T}
    
    println(io, "VoronoiMesh{$(string(Dim)),$(string(T))}")
    println(io, "   - Number of cells: $(length(msh.points))")
    println(io, "   - Number of vertices: $(length(msh.vertices))")
end


getpoint(v::VoronoiMesh) = v.points
getpoint(v::VoronoiMesh, idx) = v.points[idx]

getvertex(v::VoronoiMesh) = v.vertices
getvertex(v::VoronoiMesh, idx) = v.vertices[idx]

Base.getindex(v::VoronoiMesh) = v.cells
Base.getindex(v::VoronoiMesh, idx) = v.cells[idx]




    
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
