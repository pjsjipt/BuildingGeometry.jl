# Choping triangles with Convex polyhedrons


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
            v ∈ poly && nvin += 1
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


    # We can intersect this polygon with the triangle.
    # To do this, we need to project both the triangle and the polygon
    # on 2D plane. This is easy to do: just use the basis u⃗ and v⃗ we defined above!
end



function intersectfaceplane(f::ConvexPolygon{3,T,P}, pl{T,P,V}) where {T,P,V}
    pts = coordinates(f)

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
    for i in firstindex(s):lastindex(s)-1
        if s[i] == 0
            push!(ipts, pts[i])
        elseif s[i]*s[i+1] < 0 # There is an intersection
            # Lets calculate the intersection point
            push!(ipts, intersectpoint(pts[1], pts[2], pl))
            continue
        end
    end
    # Now we need to check the last segment
    if s[end] == 0
        push!(ipts, pts[end])
    elseif s[end]*s[begin] < 0 # AN intersection
        push!(ipts, intersectpoint(pts[end], pts[begin], pl))
    end

    return ipts
    
end

function intersectpoint(n::P, p₀::P, p₁::P, p₂::P) where{T,P<:AbstractPoint{3,T}}

    ξ = ( (p₀-p₁)⋅n ) / ( (p₂-p₁)⋅n )

    return p₁ + ξ * (p₂-p₁)
end

