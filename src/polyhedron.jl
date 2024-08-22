
import LinearAlgebra: ⋅,norm, ×
import CircularArrays: CircularVector




struct ConvexPolyhedron{T} <: GeometryPrimitive{3,T}
    "List of vertices"
    vertices::Vector{Point{3,T}}
    "Index of vertices of each face of the polyhedron"
    faces::Vector{Vector{Int}}
    vmap::Dict{Int,Int}
end

"""
`ConvexPolyhedron(vertices, faces)`

Creates a a convex polyhedron from the vertices and the faces.
The faces are specified by the indices of its vertices that make up a polygon.

The faces will be oriented so that the normal to the polygon is external to
the polyhedron.
"""
function ConvexPolyhedron(vertices::AbstractVector{Point{3,T}}, faces) where {T}
    vmap = Dict{Int,Int}()
    imap = Dict{Int,Int}()
    
    idx = 1
    for face in faces
        for i in face
            if i ∉ keys(vmap)
                vmap[i] = idx
                imap[idx] = i
                idx += 1
            end
        end
    end
    nv = length(vmap)
    verts = Vector{Point{3,T}}(undef, nv)
    for (k,v) in vmap
        verts[v] = vertices[k]
    end
    
    # Now recreate the map:
    ff = Vector{Int}[]
    for face in faces
        fi = Int[]
        for i in face
            push!(fi, vmap[i])
        end
        push!(ff, fi)
    end
    
    return ConvexPolyhedron(verts, ff, imap)
end

function Base.show(io::IO, p::ConvexPolyhedron{T}) where {T}
    
    println(io, "Polyhedron{$(string(T))}")
    println(io, "   - Number of faces: $(length(p.faces))")
    println(io, "   - Number of vertices: $(nvertices(p))")
end

vertices(p::ConvexPolyhedron) = p.vertices
GeometryBasics.coordinates(p::ConvexPolyhedron) = p.vertices

nvertices(p::ConvexPolyhedron) = length(p.vertices)
Base.length(p::ConvexPolyhedron) = length(p.vertices)

import Base.getindex
Base.getindex(p::ConvexPolyhedron, i) = ConvexPolygon(p.vertices[p.faces[i]])


nfacets(p::ConvexPolyhedron) = length(p.faces)


"""
`pnpoly(poly::ConvexPolyhedron, p::Point)`

Verify if a point is inside a polyhedron.

"""
function pnpoly(poly::ConvexPolyhedron, p::Point)

    # We will compute the outward normal of each face.
    # If the dot product between the normal and a vector
    # point from a vertex to the point is > 0, the point is
    # outside the polyhedron
    
    nf = nfacets(poly)
    for i in eachindex(poly.faces)
        pf = poly[i]
        n⃗ = normal(pf)
        u⃗ = p - vertices(pf)[begin]
        if u⃗⋅n⃗ > 0
            return false
        end
    end
    return true
end
Base.in(p::Point, poly::ConvexPolyhedron) = pnpoly(poly, p)

"""
`volume(p::ConvexPolyhedron)`

Computes the volume of a convex polyhedron.

The approach used is to find a point a inside and sum 
"""
function GeometryBasics.volume(p::ConvexPolyhedron{T}) where {T}

    # Let's choose a point: any vertex.
    pt = p.vertices[begin]
    vol = zero(T)
    for i in eachindex(p.faces)
        # Get normal and area of each face!
        face = p[i] # Convex polygon
        nrm = normalarea(face)
        A = norm(nrm)
        n⃗ = nrm ./ A
        
        # Measure the distance from point pt to the plane of the face.
        # This should be a negative number (outward normal...)
        h = (pt - vertices(face)[begin]) ⋅ n⃗
        vol -= A*h/3
    end

    return vol
                   
end

measure(p::ConvexPolyhedron) = volume(p)

"""
`centroid(p::ConvexPolyhedron)`

Compute the center of mass of a convex polyhedron
"""
function centroid(p::ConvexPolyhedron{T}) where {T}
    # As in the case of volume, we will select a point inside
    # the polyhedron (one of the vertices), split the volume in pyramids
    # that go from each face and converges to the selected point.
    # # Then we will simply add the contribution of each pyramid.
    # Let's choose a point: any vertex.
    pt = p.vertices[begin]
    xp, yp, zp = pt
    
    xc = zero(T)
    yc = zero(T)
    zc = zero(T)
    vol = zero(T)
    
    for i in eachindex(p.faces)
        # Get normal and area of each face!
        face = p[i] # Convex polygon
        x0,y0,z0 = centroid(face)
        nrm = normalarea(face)
        A = norm(nrm)
        n⃗ = nrm ./ A

        # Measure the distance from point pt to the plane of the face.
        # This should be a negative number (outward normal...)
        h = (pt - face.contour[begin]) ⋅ n⃗

        V = A*h / 3 # Volume of the pyramid
        xc += V * ( x0 + (xp-x0)/4 )
        yc += V * ( y0 + (yp-y0)/4 )
        zc += V * ( z0 + (zp-z0)/4 )

        vol += V
    end

    return Point(xc/vol, yc/vol, zc/vol)
    
end




"""
`poly2mesh(poly::ConvexPolyhedron)`

Returns a mesh of the external surface of the polyhedron. By mesh I mean
the vertices and a matrix with indices of each triangle that I decomposed the
polyheron into.

"""
function poly2mesh(poly::ConvexPolyhedron)

    nf = nfacets(poly)
    ntri = sum(length.(poly.faces)) - 2*nf
    conn = zeros(Int, ntri, 3)

    idx = 1
    for face in poly.faces
        for k in firstindex(face)+1:lastindex(face)-1
            conn[idx,1] = face[begin]
            conn[idx,2] = face[k]
            conn[idx,3] = face[k+1]
            idx += 1
        end
    end

    
    return poly.vertices, conn
                               
end

