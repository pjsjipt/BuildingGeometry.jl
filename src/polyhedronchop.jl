# Choping triangles with Convex polyhedrons
import GeometryBasics: TriangleFace

struct Plane{T,P<:AbstractPoint{3,T}, V<:AbstractPoint{3,T}}
    "Point through which the plane passes"
    p::P
    "Vector normal to the plane"
    n::V
    "Vector on the plane"
    u::V
    "Vector on the plane normal to both n and u"
    v::V
end

    

function chopwithpolyhedron(poly::ConvexPolyhedron{T}, tri::TriangleFace{P}) where {T,P<:AbstractPoint{3,T}}

    # First we will check if there is any chance of intersection
    let
        pmin,pmax = extrema(poly)
        tmin,tmax = extrema(tri)
        
        if pmin[1] > tmax[1] || pmin[2] > tmax[2] || pmin[3] > tmax[3] ||
            tmin[1] > pmax[1] || tmin[2] > pmax[2] || tmin[3] > tmax[3]
            return TriangleFace{P}[]
        end
    
        # Check if every triangle vertex is inside the Polyhedron.
        # This is another simple case
        
        nvin = 0
        for v in tri
            if v ∈ poly
                nvin += 1
            end
        end
        
        if nvin == 3 # The triangle is completely inside the polyhedron
            return TriangleFace{P}[tri]
        end
    end
    
    # Now comes the difficult part.
    # We will cut the polyhedron with the plane that contains
    # the triangle. This will result in a polygon in the plane containing the
    # triangle.
    #
    # In this plane we will intersect this polygon with the triangle
    # and then we have solved the problem
    nrm = normal(ConvexPolygon([tri[1], tri[2], tri[3]]))
    A = norm(nrm)
    n⃗ = nrm ./ A

    p₀ = tri[1]  # Point in the plane
    u⃗ = tri[2]-tri[1]
    u⃗ = u⃗ ./ norm(u⃗)
    v⃗ = n⃗ × u⃗
    plane = Plane(p₀, n⃗, u⃗, v⃗)
    
    # Let's get the polygons of each face
    segments =  [intersectfaceplane(poly[i], plane) for i in numfaces(poly)]

    # Now we need to connect the segments so that they form a polygon in the plane
    poplane = segments2polygon(segments)
    

    # We can intersect this polygon with the triangle.
    # To do this, we need to project both the triangle and the polygon
    # on 2D plane. This is easy to do: just use the basis u⃗ and v⃗ we defined above!
end



function intersectfaceplane(f::ConvexPolygon{3,T,P}, pl::Plane{T,P,V}) where {T,P,V}
    pts = coordinates(f)
    rtol = 1e-8
    p₀ = pl.p
    n⃗ = pl.n

    s = [(p-p₀)⋅n⃗ for p in pts]

    all(s .> 0) || all(s .< 0) 
    if all(s .> 0) || all(s .< 0)
        # All points of the face are above the plane or below the plane.
        # No intersection
        return P[]
    end
    ipts = P[]
    foundprev = false
    if s[begin] == 0  # Test the first point
        push!(ipts, pts[begin])
        foundprev = true
    end
    nsegs = 0 # Number of segments *in* the plane
    for i in firstindex(s):lastindex(s)-1
        if s[i+1] == 0 # The previous point has been tested. Test the end of seg.
            push!(ipts, pts[i+1])
            if foundprev
                # This means that the segment is in the plane
                nsegs += 1
            end
        elseif s[i]*s[i+1] < 0  # Test
            push!(ipts, intersectpoint(pl.n, pl.p, pts[i], pts[i+1]))
            foundprev = false
        end
    end

    # Test the last segement
    if foundprev && s[begin]==0
        nseg += 1
    end

    if nseg == 1
        # There is a single segment in the plane. We don't care about this case
        return P[]
    end
    
    if s[begin] * s[end] < 0
        push!(ipts, intersectpoint(pl.n, pl.p, pts[end], pts[begin]))
    end
        
    return ipts
    
end

"""
`intersectpoint(n::P, p₀::P, p₁::P, p₂::P)`

Given a plane defined by point `p₀` and the normal `n`, get
the intersection point with the plane of the line segment
given by points `p₁` and `p₂`.

This function assumes that there is an intersection and the segment is not
contained in the plane.
"""
function intersectpoint(n::P, p₀::P, p₁::P, p₂::P) where {T,P<:AbstractPoint{3,T}}
    Δp = p₂ - p₁
    ξ = ((p₀ - p₁)⋅n) / (Δp⋅n)
    return p₁ + ξ*Δp
end


"""
`segs2poly(segs)`

This functions takes different segments and joins them together to make a single
"""
function segments2polygon(segments::AbstractVector{Vector{P}}) where {P}

    # The intersection of a *convex* polyhedron and a plane is:
    # 1. Void
    # 2. A single point - a vertex of the polyhedron. Forget it
    # 3. A segment with 2 points
    # 4. Several points - one of the faces of the polyhedron is part of the plane
    
    # Let's start with cases 1 and 4.
    count_not_empty = 0   # Number of non empty  intersections
    count_single_int = 0  #
    for s in segments
        ns = length(s)
        if ns > 2
            return ConvexPolygon(s)
        elseif ns != 0
            count_not_empty += 1
        end
    end
    if count_not_empty == 0
        # No intersection!
        return ConvexPolygon(P[])
    end

    # The intersection 
    #
    # In the case of a segment, the seg
    
   
    
    pts = P[]
    nsegs = length(segs)
    seg_an = falses(nsegs)

    while true
        for i in eachindex(segs)
            ni = length(nsegs[i])
            
            
    for s in segs
        slen = length(s)
        if slen == 1
            
        end
    end
    
end

    
    
end
