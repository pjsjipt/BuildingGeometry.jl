
import StaticArrays: SVector
import DelimitedFiles: writedlm, readdlm
"""
`BuildingSurface(tri, points, nodes)`

A discretization of a surface, usually the outside surface of a building
undergoing wind tunnel testing. The surface is made up of a sequence if triangles,
each with a coordinate and a node.

In this context, a node is a [`NodeInfo`](@ref) object that has information on the
normal to the surface, the outside pressure tap, inside pressure tap and area.

"""
struct BuildingSurface{T,Tex,Tin}
    tri::Vector{Tri{3,T}}
    points::Vector{SVec{3,T}}
    nodes::Vector{NodeInfo{T,Tuple{Tex,Tin}}}
end

"""
`savebsurf(fname, msh)`

Save a [`BuildingSurface`](@ref) into a delimited file
"""
function savebsurf(fname, msh::BuildingSurface{T,I,I}, delim='\t') where {T,I<:Integer}

    ntri = length(msh.tri)
    tab = zeros(ntri, 13) # Vertices (9), side (2) and tag (2)
    for i in 1:ntri
        v1,v2,v3 = vertices(msh.tri[i])
        x1,y1,z1 = coordinates(v1)
        x2,y2,z2 = coordinates(v2)
        x3,y3,z3 = coordinates(v3)
        s1,s2 = nodeside(msh.nodes[i])
        t1, t2 = nodetag(msh.nodes[i])

        tab[i,:] .= [x1,y1,z1,x2,y2,z2,x3,y3,z3,s1,s2,t1,t2]
    end
    
    writedlm(fname, tab, delim)
    return
end

function loadbsurf(fname, delim='\t')
    tab = readdlm(fname, delim, Float64)
    if size(tab,2) != 13
        error("Wrong number of columns when reading BuildingSurface object. Should be 13!")
    end
    TriFace = Tri{3,Float64}
    tri = TriFace[]
    points = SVec{3,Float64}[]
    nodes = NodeInfo{3,Float64,Tuple{Int,Int}}[]

    for i in 1:size(tab,1)
        p1 = SVec(tab[i,1], tab[i,2], tab[i,3])
        p2 = SVec(tab[i,4], tab[i,5], tab[i,6])
        p3 = SVec(tab[i,7], tab[i,8], tab[i,9])

        s1 = round(Int, tab[i,10])
        s2 = round(Int, tab[i,11])
        t1 = round(Int, tab[i,12])
        t2 = round(Int, tab[i,13])

        t = TriFace(p1,p2,p3)
        p = centroid(t)
        A = area(t) * normal(t)

        node = NodeInfo(A, p, (s1, s2), (t1, t2))
        push!(tri, t)
        push!(points, p)
        push!(nodes, node)
                                 
    end
    
    return BuildingSurface(tri, points, nodes)
    
end


    

"""
`buildsurface(cad, sections; nointid=-1, tag=0)`

Discretize a section of the surface of a building into regions of influence of external and internal pressure taps. Several common situations are possible:

 1. Only external pressure taps exist
 2. Both external and internal pressure taps exist on the surface
 3. Pressure taps exist on the outside and internal pressure taps exist in some parts of the internal face of the surface.

The `cad` argument is a vector of triangles containing the geometry of the surface.
The specification of the pressure taps is given by the `sections`  argument.
This argument is a vector of named tuples.

The first element of the `sections` vector refers to the external face of the surface.
It might have 3 fields:

 1. `points`, a vector specifying the position of the external pressure taps.
 2. `tri` a vector of integers specifying the triangles in `cad` argument that make up this section of the external face of the surface.
 3. `id` An integer specifying the pressure tap id. Usually this is the index of the pressure scanner corresponding to the pressure tap.

If no internal pressure taps are present, the vector should have a single element
describing the external pressure taps. The `id` of the internal pressure tap will be
given by the keyword argument `notintid`.

Often it is useful to assign a tag to a section of the the discretization.
The keyword argument `tag` specifies the default tag used if during the section
specification, a tag is not provided.

If there are internal pressure taps, the internal face of the surface should be the
same as the external surface but it might be composed of different regions. Some of
them might have internal pressure taps and others not. The only restriction is that
the total internal face of the surface should be the same as the external surface.

Each section of the internal face of the surface has the same structure used for the
external face:

 1. `points`, a vector specifying the position of the external pressure taps.
 2. `tri` a vector of integers specifying the triangles in `cad` argument that make up this section of the external face of the surface.
 3. `id` An integer specifying the pressure tap id. Usually this is the index of the pressure scanner corresponding to the pressure tap.

Again, if `points` is not provided, no pressure taps exist.

## Example

```julia-repl
julia> epts = SVec.([(0.25, 0.25, 0), (0.75, 0.25, 0), (0.25, 0.75, 0),
                (0.75, 0.75, 0)]); # External pressure taps


julia> cad = [Tri((0,0,0),(1,0,0),(1,1,0)),
              Tri((0,0,0),(1,1,0),(0,1,0)),
              Tri((0,1,0),(1,1,0),(1,2,0)),
              Tri((0,1,0),(1,2,0),(0,2,0))];

julia> # 1. Simples case, only external pressure taps

julia> msh1 = buildsurface(cad, [(points=epts, tri=1:4, id=1:4)]; nointid=-1);

julia> # 2. Same external pressure taps and 1 internal pressure tap

julia> ipts = [SVec(0.5, 0.4, 0)]; # Internal pressure tap position

julia> msh2 = buildsurface(cad, [(points=epts, tri=1:4, id=1:4, tag=1),
                                 (points=ipts, tri=1:4, id=5, tag=2)]);

julia> # 3. Internal press. tap on the lower half and no press. tap on the upper half

julia> msh2 = buildsurface(cad, [(points=epts, tri=1:4, id=1:4, tag=1),
                                 (points=ipts, tri=1:2, id=5, tag=2),
                                 (tri=3:4, id=-1)]); # tag = 0


```
"""
function buildsurface(cad::AbstractVector{Tri{3,Float64}},
                      sections::AbstractVector; atol=atolf(Float64),nointid=-1, tag=0)

    nsecs = length(sections)

    TriFace = Tri{3,Float64}
    msh = TriFace[]
    nodes = NodeInfo{Float64,Tuple{Int,Int}}[]
    # The first section should specify the external side
    # The others, if present are internal sides
    
    # The external and internal sides should span the entire surface.
    # Different internal surfaces shouldn't overlap
    
    # Let's deal with external surface:
    esec = sections[1]
    if !haskey(esec, :points)
        # No external points. Weird...
        evor = [cad[esec.tri]]
        eid = [esec.id[1]]
        eidx = [esec.tri]
    else
        (length(esec.id) != length(esec.points)) &&
            error("`id` should have the same length as `points`")
        evor, eidx = discrsurface(cad, esec.tri, esec.points)
        eid = esec.id
    end
    if haskey(esec, :tag)
        etag = esec.tag # Tag is specified
    else
        etag = tag # Use default tag
    end
    
    if nsecs == 1  # No internal surfaces
        itag = tag  # Default tag
        for (i, m) in enumerate(evor)
            for (k,t) in enumerate(m)
                An = area(t) .* normal(t)
                tp = centroid(t)
                push!(msh, Tri(vertices(t)...))
                push!(nodes, NodeInfo(An, tp, (eid[i], nointid), (etag, itag)))
            end
        end
    else
        
        Ne = length(evor)
        
        # Let's deal with internal faces
                      itri = [collect(sections[i].tri) for i in 2:nsecs]
        etri = esec.tri
        # There should be no overlap
        nisecs = nsecs-1
        for i in 1:nisecs-1
            for k in i+1:nisecs
                if length(itri[i] ∩ itri[k]) > 0
                    error("There should be no overlap in different sections of internal surfaces!")
                end
            end
        end
        itri1 = sort(vcat(itri...))
        etri1 = sort(etri)
        
        if itri1 != etri1
            error("Internal and external surfaces should be the same")
        end
        
        # Let's start building the mesh
        for sec in sections[begin+1:end]
            if !haskey(sec, :points) # No need for voronoi. Just add the triangles
                ivor = [cad[sec.tri]]
                iid = [sec.id[1]]
                iidx = [sec.tri]
            else
                      ivor, iidx = discrsurface(cad, sec.tri, sec.points)
                iid = sec.id
            end
            if haskey(sec, :tag)
                itag = sec.tag # tag is specified
            else
                itag = tag  # Use default value of tag
            end
            
                      
            Ni = length(ivor)
            for e in 1:Ne
                for i in 1:Ni
                    intersectmesh!(msh, nodes, eid[e], evor[e], eidx[e],
                                   iid[i], ivor[i], iidx[i], etag, itag)
                end
            end
            
        end
    end

    coords = centroid.(msh)

    return BuildingSurface(msh, coords, nodes)
    
end

"""
`tri2mesh(tri)`

Create a vertices and connectivity that can be plotted using `mesh` from `Makie`.
Very simple, used for visualizing meshes. Each triangle is assumed to have independent
vertices.
"""
function tri2mesh(tri::AbstractVector{Tri{Dim,T}}) where {Dim,T}
    ntri = length(tri)
    PType = Tri{Dim,T}
    verts = zeros(T, 3, ntri*3)
    conn = zeros(Int, ntri, 3)

    k=1
    for i in 1:ntri
        vtri = vertices(tri[i])
        verts[:, k] .= vtri[1]
        verts[:, k+1] .= vtri[2]
        verts[:, k+2] .= vtri[3]
        conn[i,:] .= (k,k+1,k+2)
        k += 3
    end
    return verts, conn
end
tri2mesh(m::BuildingSurface) = tri2mesh(m.tri)

"""
`floor_mesh(tri, idx, msh)`

Create a [`BuildingSurface`](@ref) for each floor of a building.

## Arguments
 * `tri` Triangles that make up the surface of the floor
 * `idx` Index of each `tri` corresponding to the global mesh
 * `msh` Global Building discretization

The building is initially discretized into a global [`BuildingSurface`](@ref)
structure given by argument `msh`. Each floor (or slice) is made up of triangles
obtained by the method [`slicemesh`](@ref).

"""
function floor_mesh(tri, idx, msh)
    A = [area(t) .* normal(t) for t in tri]
    coords = centroid.(tri)
    iex = [n.side[1] for n in msh.nodes[idx]]
    iin = [n.side[2] for n in msh.nodes[idx]]
    etag = nodetag.(msh.nodes, 1)
    itag = nodetag.(msh.nodes, 1)
    
    nodes = [NodeInfo(A[i], coords[i], (iex[i], iin[i]), (etag[i],itag[i]))
             for i in eachindex(A)]
    
    return BuildingSurface(tri, coords, nodes)
    
end

"""
`buildingslice(msh, p)`
`buildingslice(msh, z)`
`buildingslice(msh, nslices, pa, pb)`

Slice the global mesh of a building into slices. The slices can be specified in
different ways depending on the arguments types.

## Arguments

 * `msh` Global mesh of a building
 * `p` A vector of points. Each slice i is specified by the region `p[i]` to `p[i+1]`
 * `z` The slices can also be specified by the coordinates in the z direction
 * `nslices`, `pa`, `pb`. The building is sliced into equal `nslices` from point `pa` to point `pb`.
"""
function buildingslice(msh::BuildingSurface{T},
                    p::AbstractVector{SVec{3,T}};  atol=atolf(Float64)) where {T}
    trilst, triidx = slicemesh(msh.tri, p; atol=atol)
    
    return [floor_mesh(trilst[i], triidx[i], msh) for i in eachindex(trilst)]
end


function buildingslice(msh::BuildingSurface{T}, z::AbstractVector{T};
                       x=0, y=0, atol=atolf(Float64)) where {T}
    
    trilst, triidx = slicemesh(msh.tri, z; x=x, y=y, atol=atol)
    return [floor_mesh(trilst[i], triidx[i], msh) for i in eachindex(trilst)]
end

    
function buildingslice(msh::BuildingSurface{T}, nslices::Integer,
                   pa::SVec{3,T}, pb::SVec{3,T};   atol=atolf(Float64)) where {T}

    trilst, triidx = slicemesh(msh.tri, pa, pb; atol=atol)
    
    return [floor_mesh(trilst[i], triidx[i], msh) for i in eachindex(trilst)]
end


"""
`mergemeshes(mshlst)`

Each surface is discretized independently resulting in a [`BuildingSurface`](@ref)
object. This contains information on internal and external pressure taps, triangles
and positions. This function merges the meshes into a single [`BuildingSurface`](@ref)
that can be used to view results in 3d plots or slice the building into floors.
"""
function mergemeshes(mshlst::AbstractVector{<:BuildingSurface{T}}) where {T}
    tri = copy(mshlst[begin].tri)
    points = copy(mshlst[begin].points)
    nodes = copy(mshlst[begin].nodes)
    for msh in mshlst[begin+1:end]
        append!(tri, msh.tri)
        append!(points, msh.points)
        append!(nodes, msh.nodes)
    end
    
    return BuildingSurface(tri, points, nodes)
    
end

mergemeshes(mshlst...) = mergemeshes(collect(mshlst))

