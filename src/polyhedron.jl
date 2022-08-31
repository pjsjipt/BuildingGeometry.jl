
import LinearAlgebra: ⋅,norm, ×

import GeometryBasics
import GeometryBasics: Polytope
import GeometryBasics: TriangleFace

export ConvexPolyhedron, numfaces, poly2mesh, volume, numvertices, Rect



struct ConvexPolyhedron{T,P<:AbstractPoint{3,T},VL<:AbstractVector{P}}
    "List of vertices"
    vertices::VL
    "Index of each individual vertex"
    vlist::Vector{Int}
    "Index of vertices of each face of the polyhedron"
    faces::Vector{Vector{Int}}
end

function ConvexPolyhedron(vertices, faces)
    vlist = Set{Int}()
    
    nf = length(faces)
    ff = Vector{Int}[]
    for face in faces
        fi = Int[]
        for i in face
            push!(vlist, i)
            push!(fi,i)
        end
        push!(ff, fi)
    end

    return ConvexPolyhedron(vertices, collect(vlist), faces) 
end

function Base.show(io::IO, p::ConvexPolyhedron{T}) where {T}
    
    println(io, "Polyhedron{$(string(T))}")
    println(io, "   - Number of faces: $(length(p.faces))")
    println(io, "   - Number of vertices: $(length(p.vlist))")
end


import Base.getindex
Base.getindex(p::ConvexPolyhedron, i) = ConvexPolygon(p.vertices[p.faces[i]])

getvertex(p::ConvexPolyhedron) = p.vertices[p.vlist]
getvertex(p::ConvexPolyhedron, i) = p.vertices[p.vlist[i]]

numfaces(p::ConvexPolyhedron) = length(p.faces)
numvertices(p::ConvexPolyhedron) = length(p.vlist)

import GeometryBasics.Rect

Rect(p::ConvexPolyhedron) = Rect(p.vertices[p.vlist])

"""
`pnpoly(poly::ConvexPolyhedron, p::Point)`

Verify if a point is inside a polyhedron.

"""
function pnpoly(poly::ConvexPolyhedron, p::Point)

    # We will compute the outward normal of each face.
    # If the dot product between the normal and a vector
    # point from a vertex to the point is > 0, the point is
    # outside the polyhedron
    
    nf = numfaces(poly)
    for i in eachindex(poly.faces)
        pf = poly[i]
        n⃗ = normal(pf)
        u⃗ = p - poly.vertices[poly.faces[i][begin]]
        if u⃗⋅n⃗ > 0
            return false
        end
    end
    return true
end
Base.in(p::Point, poly::ConvexPolyhedron) = pnpoly(poly, p)

import GeometryBasics: volume
"""
`volume(p::ConvexPolyhedron)`

Computes the volume of a convex polyhedron.

The approach used is to find a point a inside and sum 
"""
function volume(p::ConvexPolyhedron{T}) where {T}

    # Let's choose a point: any vertex.
    pt = p.vertices[p.vlist[begin]]
    vol = zero(T)
    for i in eachindex(p.faces)
        # Get normal and area of each face!
        face = p[i] # Convex polygon
        nrm = normal(face)
        A = norm(nrm)
        n⃗ = nrm ./ A

        # Measure the distance from point pt to the plane of the face.
        # This should be a negative number (outward normal...)
        h = (pt - face.contour[begin]) ⋅ n⃗
        vol -= A*h/3
    end

    return vol
                   
end

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
    xp, yp, zp = pt = p.vertices[p.vlist[begin]]

    xc = zero(T)
    yc = zero(T)
    zc = zero(T)
    vol = zero(T)
    
    for i in eachindex(p.faces)
        # Get normal and area of each face!
        face = p[i] # Convex polygon
        x0,y0,z0 = centroid(face)
        nrm = normal(face)
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

import Base.extrema
function extrema(p::ConvexPolyhedron)
    nv = length(p.vlist)
    iv = p.vlist
    v = p.vertices
    p1 = v[iv[begin]]

    xmax = xmin = p1[1]
    ymax = ymin = p1[2]
    zmax = zmin = p1[3]

    for i in iv
        x,y,z = v[i]
        if x > xmax
            xmax = x
        elseif x < xmin
            xmin = x
        end
        if y > ymax
            ymax = y
        elseif y < ymin
            ymin = y
        end
        if z > zmax
            zmax = z
        elseif z < zmin
            zmin = z
        end
        
    end
    return Point3(xmin, ymin, zmin), Point3(xmax, ymax, zmax)
end

function extrema(p::TriangleFace)

    xmin = min(p[1][1], p[2][1], p[3][1])
    xmax = max(p[1][1], p[2][1], p[3][1])
    ymin = min(p[1][2], p[2][2], p[3][2])
    ymax = max(p[1][2], p[2][2], p[3][2])
    zmin = min(p[1][3], p[2][3], p[3][3])
    zmax = max(p[1][3], p[2][3], p[3][3])
    return Point3(xmin, ymin, zmin), Point3(xmax, ymax, zmax)
end
    
normal(tri::TriangleFace) = (tri[3]-tri[1]) × (tri[2]-tri[1]) ./ 2
area(tri::TriangleFace) = norm(normal(tri))


"""
`poly2mesh(poly::ConvexPolyhedron)`

Returns a mesh of the external surface of the polyhedron. By mesh I mean
the vertices and a matrix with indices of each triangle that I decomposed the
polyheron into.

"""
function poly2mesh(poly::ConvexPolyhedron)

    nf = numfaces(poly)
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
