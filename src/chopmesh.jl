# Choping triangles with Convex polyhedrons
import StaticArrays: SVector


import StaticArrays: SVector

function plane_side(p, pl_p0, pl_n)

    s = (p - pl_p0) ⋅ pl_n
    if isapproxzero(s)
        return 0
    elseif s > zero(s)
        return 1
    else
        return -1
    end
end


plane_side(p, plane::Plane) = plane_side(p, plane.p, plane.n)
@inline idxn(i) = (i % 3) + 1
@inline idxp(i) = (i+1) % 3 + 1
function cut_with_plane(tri::TriangleFace, plane::Plane)

    out = typeof(tri)[]
    p = coordinates(tri)
    TT = typeof(tri)
    p0 = plane.p
    n = normal(plane)
    
    s = plane_side.(p, Ref(p0), Ref(n))

    # A plane can be defined by a point and its normal
    # Funcion `plane_side` can be used to check if a point is
    # 1. On the side pointed by the normal (+)
    # 2. On the side opposite to the normal (-)
    # 3. On the plane itself.
    # If all points are either on the plane or on the + side,
    # then nothing should be returned. If all points are on
    # the - side or on the plane, the triangle itself should be returned.
    
    

    if all(s .>= 0)
        # Every vertex of the triangle is on the + side of the triangle. Nothing!
        return out
    elseif all(s .<= 0)
        # Every vertex is on the - side of the plane: return the triangle itself
        push!(out, tri)
        return out
    end

    
    # Now we have a more complex case. Let's check edge by edge
    # There are 3 possible situations
    # 1. Two vertices on the + side and 1 in - side
    # 2. Two vertices on the - side and 1 in + side
    # 3. One vertex in the + side, one in the - side and one on the plane

    s1 = (s[1], s[2], s[3], s[1])

    if any(s .== 0)
        idx = findfirst(isequal(0), s)
        inx = idxn(idx)
        ipr = idxp(idx)
        pi = intersectpoint(n, p0, p[inx], p[ipr])
        if s[inx] < 0
            push!(out, TT(p[idx], p[inx], pi))
        else
            push!(out, TT(p[idx], pi, p[inx]))
        end
           
    elseif sum(s .== 1) == 2
        idx = findfirst(isequal(-1), s)
        inx = idxn(idx)
        ipr = idxp(idx)
        pi1 = intersectpoint(n, p0, p[idx], p[inx])
        pi2 = intersectpoint(n, p0, p[idx], p[ipr])

        push!(out, TT(p[idx], pi1, pi2))
    else
        idx = findfirst(isequal(1), s)
        inx = idxn(idx)
        ipr = idxp(idx)
        pi1 = intersectpoint(n, p0, p[idx], p[inx])
        pi2 = intersectpoint(n, p0, p[idx], p[ipr])

        push!(out, TT(pi1, p[inx], pi2))
        push!(out, TT(pi2, p[inx], p[ipr]))
    end

    return out
    
end
"""
`cut_with_plane(pts, p0, n, circ; atol=1e-8)`

Given a polygon or set of linked lines characterized by a vector of points,
this function will cut it using a plane characterize by a point `p0` and normal `n`.

The function will return the points inside the plane. If the plane crosses the
polygons or lines, new segments that consist of the intersection of the polygon
and the plane will be added.

This function is used throughout this package for many types of operations.

"""
function cut_with_plane(pts::Union{AbstractVector{Point{3,T}},NTuple{N,Point{3,T}}},
                        p0, n, circ=true; atol=1e-8) where{T,N}

    npts = length(pts)
    out = Point{3,T}[]
    if length(pts) < 1
        return out
    end
    zz = zero(T)
    p1 = pts[begin]
    slast = (p1 - p0) ⋅ n
    if isapprox(slast, zz, atol=atol)
        slast = zz
    end
    if slast <= 0
        # The first point is inside the plane. Point should be conserved!
        if !circ
            push!(out, p1)
        end
    end

    if circ
        idx = [2:npts; 1]
    else
        idx = collect(2:npts)
    end
 
    for i in idx
        p2 = pts[i]
        snew = (p2 - p0) ⋅ n
        if isapprox(snew, zz, atol=atol)
            snew = zz
        end
       
        if snew == 0.0
            push!(out, p2)
        elseif snew * slast  < 0
            px = intersectpoint(n, p0, p1, p2)
            push!(out, px)
            if snew < 0
                push!(out, p2)
            end
        elseif snew < 0
            push!(out, p2)
        end
        p1 = p2
        slast = snew
    end
    return out
end

"""
`chopwithpolyhedron(poly, tri; atol=1e-8)`

Uses a polyhedron to chop a triangle, meaning that the function
will return parts of the triangle that are insided the polygon.

First, this function checks if the boundaing boxes of polyhedron and triangle
overlap. If they don't, the result is empty.

If every vertex of the triangle is inside the polyhedron, just return the triangle.

Then, for each face of the polyhedron, the function calls [`cut_with_plane`](@ref)
chopping away the parts of the triangle that are outside the polyhedron.

The function returns a list of triangles that are inside the polyhedron.

"""
function chopwithpolyhedron(poly::ConvexPolyhedron{T},
                            tri::TriangleFace) where {T}
    
    # First we will check if there is any chance of intersection
    TT = typeof(tri)
    let
        pbox = boundingbox(poly)
        tbox = boundingbox(tri)
        if pbox.min > tbox.max || pbox.min > tbox.max || pbox.min > tbox.max ||
            tbox.min > pbox.max || tbox.min > pbox.max || tbox.min > tbox.max
            return TT[]
        end
        # Check if every triangle vertex is inside the Polyhedron.
        # This is another simple case

        nvin = 0
        for v in coordinates(tri)
            if v ∈ poly
                nvin += 1
            end
        end
        
        if nvin == 3 # The triangle is completely inside the polyhedron
            return [TT(coordinates(tri)...)]
        end
    end
    # Here is the tricky part.
    # Each face of the *convex* polyhedron is part of a plane that divides the
    # entire space into two parts. Outside and inside. Points that are insider
    # the polyhedron is inside every polygonal face.
    #
    # So to chop the triangle using the polyhedron, we will use the same idea:
    # 1. We will start with a polygon with 3 vertices - the triangle.
    # 2. We will sweep this polygon. If a vertex is outside the plane,
    #    this vertex we need to see if the previous vert

    
    nf = nfacets(poly)
    t1 = [tri]
    for i in 1:nf
        face = poly[i]
        pl = Plane(face) 
        t2 = TT[]
        for t in t1
            tx = cut_with_plane(t, pl)
            for x in tx
                push!(t2, x)
            end
        end
        t1 = t2
    end
    return t1
end

    

"""
`intersectpoint(n::V, p₀::P, p₁::P, p₂::P)`

Given a plane defined by point `p₀` and the normal `n`, get
the intersection point with the plane of the line segment
given by points `p₁` and `p₂`.

This function assumes that there is an intersection and the segment is not
contained in the plane.
"""
function intersectpoint(n::Vec, p₀::Point, p₁::Point, p₂::Point)
    Δp = p₂ - p₁
    ξ = ((p₀ - p₁)⋅n) / (Δp⋅n)
    return p₁ + ξ*Δp
end



    

