# Choping triangles with Convex polyhedrons


export cut_with_plane, chopwithpolyhedron

function cut_with_plane(pts::Vector{P}, p0, n, circ=true; atol=1e-8) where{T,P<:AbstractPoint{T}}

    
    npts = length(pts)
    out = P[]
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


function chopwithpolyhedron(poly::ConvexPolyhedron{T}, tri::TriangleFace{P}; atol=1e-8) where {T,P<:AbstractPoint{3,T}}

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

    # Here is the tricky part.
    # Each face of the *convex* polyhedron is part of a plane that divides the
    # entire space into two parts. Outside and inside. Points that are insider
    # the polyhedron is inside every polygonal face.
    #
    # So to chop the triangle using the polyhedron, we will use the same idea:
    # 1. We will start with a polygon with 3 vertices - the triangle.
    # 2. We will sweep this polygon. If a vertex is outside the plane,
    #    this vertex we need to see if the previous vert

    pts = [tri[1], tri[2], tri[3]]
    
    nf = numfaces(poly)
    
    for i in 1:nf
        face = poly[i]
        n = normal(face)
        p0 = coordinates(face)[1]
        pts = cut_with_plane(pts, p0, n, atol=atol)
    end
    return TriangleFace{P}[TriangleFace(pts[1], pts[i], pts[i+1]) for i in 2:length(pts)-1]
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



    

