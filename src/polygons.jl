

using StaticArrays

struct ConvexPolygon{M<:Manifold,C<:CRS,R<:Ring{M,C}} <: Polygon{M,C}
    contour::R
end
"""
`ConvexPolygon(contour)`
`ConvexPolygon(x,y)`
`ConvexPolygon(x,y,z)`
`ConvexPolygon(pts)`

Create a convex polygon. 


"""

ConvexPolygon(x::AbstractVector, y::AbstractVector) =
    ConvexPolygon([Point(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    ConvexPolygon([Point(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

function ConvexPolygon(pts::AbstractVector{P}) where {P<:Point}
    ConvexPolygon(Ring(pts))
end

ConvexPolygon(pts::AbstractVector{TP}) where {TP<:Tuple} = ConvexPolygon(Point.(pts))
ConvexPolygon(pts::Vararg{P}) where {P<:Point} = ConvexPolygon(collect(pts))

Meshes.vertices(p::ConvexPolygon) = vertices(p.contour)
Meshes.nvertices(p::ConvexPolygon) = nvertices(p.contour)
Meshes.rings(p::ConvexPolygon) = [p.contour]
                                     
"""
`normalarea(p)`

Compute the area vector of a convex polygon. By area vector, it is meant the outer
 vector whose length is the area of the convex polygon.
"""
function normalarea(p::ConvexPolygon)
    pts = vertices(p)
    sum( (pts[i-1]-pts[begin]) × (pts[i]-pts[begin])
         for i in firstindex(pts)+1:lastindex(pts) ) ./ 2
end

function normal_(p::ConvexPolygon)
    na = normalarea(p)
    return normalize(na)
end


#function Meshes.normal(p::ConvexPolygon{2})
#    n = normal_(p)
#    return n / abs(n)
#end


Meshes.area(p::ConvexPolygon) = norm(normalarea(p))
Meshes.measure(p::ConvexPolygon) = area(p)

function Meshes.normal(p::ConvexPolygon)
    pts = vertices(p)
    unormalize(  sum( ucross(pts[i-1]-pts[begin], pts[i]-pts[begin])
                      for i in firstindex(pts)+1:lastindex(pts) )    )
    
end

"""
`pnpoly(poly, p)`

Determine if a point is inside a polygon. If a point is on the boundary of the polygon
this function might return `true` or `false` depending on floating point errors and
such.
"""
function pnpoly(poly::ConvexPolygon{𝔼{2}}, p::Point)

    v = vertices(poly)

    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = ustrip(udot( (v₂ - v₁),  (p - v₁) ))
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end

function pnpoly(poly::ConvexPolygon{𝔼{3}}, p::Point)

    
    v = vertices(poly)
    n = normal(poly)

    pl = Plane(v[1], n) # Plane

    # Check whether the point is in the plane of polygon
    if p ∉ pl
        return false
    end

    
    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = ustrip(udot( (v₂ - v₁),  (p - v₁) ))
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end
                
Meshes.in(pt::Point, p::ConvexPolygon) = pnpoly(p, pt)


function Meshes.centroid(p::ConvexPolygon)
    pts = vertices(p)

    p1 = to(pts[begin])  # Let's work with vectors
    up1 = ustrip.(p1)
    un = unit(p1[1])
    T = typeof(ustrip(p1[1])) # Type (float)
    N = length(p1)
    A = zero(T)
    CA = zero(ustrip.(p1))
    for i in firstindex(pts)+1:lastindex(pts)-1
        u = ustrip.(pts[i]-pts[1])
        v = ustrip.(pts[i+1]-pts[1])
        
        Ai = norm(u × v) / 2  # Area of the triangle
        A += Ai
        # Compute the centroid of the triangle: mean of the coordinates
        Ci = ustrip.( p1 + to(pts[i]) + to(pts[i+1]) ) ./ 3
        CA += Ai .* Ci
    end
    return Point(  ((CA./A) .* un )... ) 
            
end



function Meshes.simplexify(poly::ConvexPolygon)

    pts = vertices(poly)

    conn = [connect((1,i,i+1)) for i in firstindex(pts)+1:lastindex(pts)-1]
    return SimpleMesh(pts, conn)
  
    
end

Meshes.discretize(poly::ConvexPolygon, method=nothing) = simplexify(poly)

plane(p::ConvexPolygon) = Plane(vertices(p)[begin], normal(p))




"""
`poly2mesh(poly)`

Creates a triangle mesh from a convex polygon. 
"""
function poly2mesh(poly::ConvexPolygon)
    pts = vertices(poly)
    
   
    ntri = length(pts) - 2
    conn = zeros(Int, ntri, 3)

    idx = 1
    for i in firstindex(pts)+1:lastindex(pts)-1
        conn[idx,1] = firstindex(pts)
        conn[idx,2] = i
        conn[idx,3] = i+1
        idx += 1
    end

    return pts, conn
end


