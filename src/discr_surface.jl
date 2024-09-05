export discrsurface
# Surface discretization
import Unitful: m
m² = m*m
"""
`discrsurface(tri, idx, pts)`
`discrsurface(tri, pts)`

Given a surface made up of triangles in vector `tri` this function will
find the influence region of each point in vector `pts` assumed to be located
on the surface. If only parts the triangle mesh make up the surface, parameter
`idx` can be used to select the triangles that should be considered.

The function computes the 3D Voronoi diagram of each point. Then uses the resulting
polyhedrons to chop the triangle mesh into regions of influence of each point.

Region of influence here is the parts of the triangle mesh that is closest to a point
in `pts` than anyother point in `pts`.

Care should be taken to choose the points well: a region on another side of the surface
could potentially be closest to a point if the points are not well distributed.

The function returns the triangles that make up each region of influence. It returns, as well, the index of the original triangle from where it was created.
"""
function discrsurface(tri, idx::AbstractVector{<:Integer},
                      pts::AbstractVector{<:Point};
                      bbox=nothing, nd=8, return_vor=false)

    TriFace = eltype(tri)
    
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

    vor = voronoi3d(pts, bbox=bbox)

    if return_vor
        return vor 
    end

    
    # Chop each triangle with every polyhedron.
    # Each node of `pts` corresponds to a volume (polyhedron). We will
    # Get every triangle on the surface mesh and 
    trivor = Vector{TriFace}[]
    tidx = Vector{Int}[]
    for vol in vor.cells
        id = Int[]
        trim = TriFace[]
        for i in idx
            
            t = tri[i]
            m = chopwithpolyhedron(vol, t)
            nm = length(m)
            if nm > 0
                for ti in m
                    aa = area(ti)
                    if ispositive(aa)
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

discrsurface(tri, pts::AbstractVector{<:Point};
             bbox=nothing, nd=8, return_vor=false) =
                 discrsurface(tri, 1:length(tri), pts;
                              bbox=bbox,nd=nd, return_vor=return_vor)



"""
`slicemesh(m, p)`
`slicemesh(m, z; x=0, y=0)`
`slicemesh(m, nslices, pa, pb;  rtol=1e-8)`

Slice a mesh in slices defined by a vector of points `p`. Since in most cases
the slices will be along the z axis, there is a second method that takes that into
account. There is a method for slicing `nslices` from `pa` to `pb`.

## Arguments

 * `m`: Vector of geometric shapes. Should be a convex polygon or something close
 * `p`: Vector of points defining the planes used for slicing the mesh
 * `z`: Alternatively, the height at which the mesh should be sliced
 * `rtol`: relative tolerance admitted.

"""
function slicemesh(m::AbstractVector{<:Triangle},
                   p::AbstractVector{<:Point})
    TriFace = eltype(m)
    mshlst = Vector{TriFace}[]
    mshidx = Vector{Int}[]
    # We will analyze each slice
    for i in firstindex(p)+1:lastindex(p)
        idx = Int[] # Id's of each geometry
        mshi = TriFace[] # We will decompose this into triangles
        p₁ = p[i-1]
        p₂ = p[i]
        n⃗₂ = Meshes.unormalize(p₂ - p₁) 
        n⃗₁ = -n⃗₂
        pl1 = Plane(p₁, n⃗₁) # Plane 1
        pl2 = Plane(p₂, n⃗₂) # Plane 1
        
        for (k,f) in enumerate(m)
            v = vertices(f)
            # Check if there are vertices or edges in the slice in question
            
            all(isnegative(udot(u-p₁,n⃗₂)) for u in v) && continue # Vertices below. Next!
            all(ispositive(udot(u-p₂, n⃗₂)) for u in v) && continue # Vertices above. Next!
            
            # We have an intersection
            # Chop with plane 1
            tri1 = cut_with_plane(f, pl1)
            # Chop with plane 2. Now we might have more than 1 triangle
            tri2 = TriFace[]
            for tri in tri1
                xtri = cut_with_plane(tri, pl2)
                for t in xtri
                    push!(tri2, t)
                end
            end
            
            # Let's store this information
            for t in tri2
                push!(mshi, t)
                push!(idx, k)
            end
            
        end
        push!(mshlst, mshi)
        push!(mshidx, idx)
    end

    return mshlst, mshidx
end

function slicemesh(msh::AbstractVector{<:Triangle}, z::AbstractVector;
                   x=0m, y=0m)
    p = Point.(x, y, z)
    return slicemesh(msh, p)
end


function slicemesh(m::AbstractVector{<:Triangle}, nslices::Integer,
                   pa::Point, pb::Point)
    ξ = range(0.0, 1.0, length=nslices+1)
    u⃗ = pb-pa
    p = [pa + ξᵢ * u⃗ for ξᵢ in ξ]
    return slicemesh(m, p)
end

