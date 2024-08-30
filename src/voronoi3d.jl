
export voronoi3d,  makebbox, VoronoiMesh
export getpoint, getvertex, npoints


function bboxaux!(xlims, p)
    x = repeat([xlims[1], xlims[2]], inner=p[1], outer=p[2])
    p[1] = p[1] * 2
    p[2] = p[2] ÷ 2
    return x
end


function makebbox(ranges...)
    nr = length(ranges)
    p = [1,2^(nr-1)] # Outer and inner initial
    return map(x->bboxaux!(x, p), ranges)
    
end

"""
`voronoi(x,y,z; bbox,nd,auxpts)`
`voronoi(pts; bbox,nd,auxpts)`

Calls `scipy.spatial.Voronoi` to compute the 3D Voronoi diagram of a collection
of 3D points.

Before computing the Voronoi diagram, auxiliary points, are added and the corresponding
volumes are excluded. `scipy.spatial.Voronoi` call `QHull` library and the points
at infinity are all the same which is a difficulty. If a bounding box is provieded,
the extremities are auxilary points which help creating the diagram without the infinity
points. If no bounding box is provided, parameter `nd` is used to create it using `nd` × the largest separation of the points. Further auxiliary points can be provided.

The function returns a [`VoronoiMesh`](@ref) object which contains every polyhedron
of the Voronoi diagram and the connectivity - between nodes and vertices.
"""
function voronoi3d(x,y,z; bbox=nothing, nd=8, auxpts=nothing)
    length(x) == length(y) == length(z) || error("Vectors x, y and z should have the same length")
    
    if isnothing(bbox)
        bbox = let
            xmin,xmax = extrema(x)
            ymin,ymax = extrema(y)
            zmin,zmax = extrema(z)
            
            Δ = max(xmax-xmin, ymax-ymin, zmax-zmin)  # A length for the domain
            # Let's make a big bounding box
            Box( Point(xmin-nd*Δ, ymin-nd*Δ, zmin-nd*Δ),
                 Point(xmax+nd*Δ, ymax+nd*Δ, zmax+nd*Δ) )
        end
        
    end
    
    #=
    The coordinates of the bounding box will make up 4 extra points
    at the beggining of the point list. This is done so that the
    volumes that contain the infinity point can be excluded.
    =#
    bx, by, bz  = let
        (xmin,ymin,zmin), (xmax,ymax,zmax) = bbox.min, bbox.max
        bx1, by1, bz1 = makebbox((xmin,xmax),
                                 (ymin,ymax),
                                 (zmin,zmax))
        
        xm = (xmin+xmax)/2
        ym = (ymin+ymax)/2
        zm = (zmin+zmax)/2
        bx2 = [xmin, xmax, xm, xm, xm, xm]
        by2 = [ym, ym, ymin, ymax, ym, ym]
        bz2 = [zm, zm, zm, zm, zmin, zmax]
        if !isnothing(auxpts)
            [bx1; bx2; auxpts.x], [by1; by2; auxpts.y], [bz1; bz2; auxpts.z]
        else
            [bx1; bx2], [by1; by2], [bz1; bz2]
        end

    end
    nbb = length(bx)

    px = [bx; x]
    py = [by; y]
    pz = [bz; z]
    
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
        pface = ConvexPolygon(CircularVector(v[vi]))
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

    # Build Vertex connectivity:
    nverts = length(v)
    vconn = [Int[] for i in 1:nverts]
    vpconn = [Int[] for i in 1:nverts]
    for (i, f) in enumerate(plst)
        # i is the point. f is a vector with vertices
        for fi in f
            nf = length(fi)
            for k in 1:nf
                vk = fi[k]
                for j in k+1:nf
                    vj = fi[k]
                    if vj ∉ vconn[vk]  # We have already considered this pair
                        push!(vconn[vk], vj)
                    end
                    if vk ∉ vconn[vj]
                        push!(vconn[vj], vk)
                    end
                end
                if i ∉ vpconn[vk]
                    push!(vpconn[vk], i)
                end
            end
        end
    end

    return VoronoiMesh(p, v, cells, conn, vconn, vpconn)
    
end

function voronoi3d(pts::AbstractVector{<:Point{3}}; bbox=nothing, nd=8, auxpts=nothing)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    z = [p[3] for p in pts]

    if !isnothing(auxpts)
        xa = [p[1] for p in auxpts]
        ya = [p[2] for p in auxpts]
        za = [p[3] for p in auxpts]
        auxpts = (x=xa, y=ya, z=za)
    end
    
    return voronoi3d(x, y, z; bbox=bbox, nd=nd, auxpts=auxpts)
end

             
    
"""
`VoronoiMesh(p,v,conn,vconn,vpconn)`

Stores information related to the Voronoi diagram of a set of points.
This mesh stores the points used to generate the diagram, the vertices
of the diagram, the convex polyhedrons of each volume of the Voronoi diagram
and connectivity. It has connectivity information related to

 * Points Neighboring a point
 * Vertices neighboring a vertex
 * Points neighboring a vertex
"""
struct VoronoiMesh{Dim,T,VP<:AbstractVector{Point{Dim,T}},
                   VV<:AbstractVector{Point{Dim,T}},CP,
                   VPO<:AbstractVector{CP}}
    "Points used to calculate the Voronoi diagram"
    points::VP
    "Vertices of the Voronoi diagram"
    vertices::VV
    "Individual cells of the Voronoi diagram"
    cells::VPO
    "Point connectivity"
    conn::Vector{Vector{Int}}
    "Vertex connectivity"
    vconn::Vector{Vector{Int}}
    "Vertex point connectivity"
    vpconn::Vector{Vector{Int}}
end

function Base.show(io::IO, msh::VoronoiMesh{Dim,T}) where {Dim,T}
    
    println(io, "VoronoiMesh{$(string(Dim)),$(string(T))}")
    println(io, "   - Number of cells: $(length(msh.points))")
    println(io, "   - Number of vertices: $(length(msh.vertices))")
end

GeometryBasics.coordinates(v::VoronoiMesh) = v.vertices
Base.length(v::VoronoiMesh) = length(v.vertices)
vertices(v::VoronoiMesh) = v.vertices
nvertices(v::VoronoiMesh) = length(v.vertices)


npoints(v::VoronoiMesh) = length(points)

points(v::VoronoiMesh) = v.points


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
