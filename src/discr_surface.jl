export discrsurface
# Surface discretization

function discrsurface(tri, idx::AbstractVector{<:Integer},
                      pts::AbstractVector{Point{3,Float64}};
                      rtol=1e-8, bbox=nothing, nd=8)

    TriFace = Triangle{3,Float64,SVector{3,Point{3,Float64}}}
    
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
    trivor = Vector{TriFace}[]
    tidx = Vector{Int}[]
    for vol in vor.cells
        id = Int[]
        trim = TriFace[]
        for i in idx
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

discrsurface(tri, pts::AbstractVector{Point{3,Float64}};
             rtol=1e-8, bbox=nothing, nd=8) = discrsurface(tri, 1:length(tri), pts;
                                                           rtol=rtol, bbox=bbox,nd=nd)



"""
`slicemesh(m, p; rtol=1e-8)`
`slicemesh(m, z; x=0, y=0, rtol=1e-8)`
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
function slicemesh(m::AbstractVector{P},
                   p::AbstractVector{Point{3,T}}; rtol=1e-8) where {P,T}
    mshlst = Vector{Triangle{3,T}}[]
    mshidx = Vector{Int}[]
    # We will analyze each slice
    for i in firstindex(p)+1:lastindex(p)
        idx = Int[] # Id's of each geometry
        mshi = Triangle{3,T}[] # We will decompose this into triangles
        p₁ = p[i-1]
        p₂ = p[i]
        L = norm(p₂-p₁)
        n⃗₂ = (p₂ - p₁) ./ L
        n⃗₁ = -n⃗₂
        atol = rtol*L
        for (k,f) in enumerate(m)
            v = vertices(f)
            # Check if there are vertices or edges in the slice in question
            
            all((u-p₁) ⋅ n⃗₂ < 0 for u in v) && continue # Vertices below. Next!
            all((u-p₂) ⋅ n⃗₂ > 0 for u in v) && continue # Vertices above. Next!

            # We have an intersection
            pts = cut_with_plane(v, p₁, n⃗₁, true; atol=atol)
            pts = cut_with_plane(pts, p₂, n⃗₂, true; atol=atol)

            # Generate triangles:
            npts = length(pts)
            for i in 2:npts-1
                push!(mshi, Triangle(pts[1], pts[i], pts[i+1]))
                push!(idx, k)
            end
        end
        push!(mshlst, mshi)
        push!(mshidx, idx)
    end

    return mshlst, mshidx
end

function slicemesh(m::AbstractVector{P}, z::AbstractVector{T};
                   x=0, y=0, rtol=1e-8) where {P,T}
    p = Point{3,T}.(x, y, z)
    return slicemesh(m, p; rtol=rtol)
end


function slicemesh(m::AbstractVector{P}, nslices::Integer,
                   pa::Point{3,T}, pb::Point{3,T};  rtol=1e-8) where {P,T}
    ξ = range(0.0, 1.0, length=nslices+1)
    u⃗ = pb-pa
    p = [pa + ξᵢ * u⃗ for ξᵢ in ξ]
    return slicemesh(m, p; rtol=rtol)
end

