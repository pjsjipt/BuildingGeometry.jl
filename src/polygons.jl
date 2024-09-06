

using StaticArrays
import CircularArrays: CircularVector
struct ConvexPolygon{Dim,T} <: AbstractBuildGeom
    contour::CircularVector{SVec{Dim,T}}
end
"""
`ConvexPolygon(contour)`
`ConvexPolygon(x,y)`
`ConvexPolygon(x,y,z)`
`ConvexPolygon(pts)`

Create a convex polygon. 


"""

ConvexPolygon(x::AbstractVector, y::AbstractVector) =
    ConvexPolygon([SVec(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    ConvexPolygon([SVec(xx, yy, zz) for (xx,yy,zz) in zip(x,y,z)])

function ConvexPolygon(pts::AbstractVector{P}) where {P<:SVec}
    ConvexPolygon(CircularVector(pts))
end

ConvexPolygon(pts::AbstractVector{TP}) where {TP<:Tuple} = ConvexPolygon(SVec.(pts))
ConvexPolygon(pts::Vararg{P}) where {P<:SVec} = ConvexPolygon(collect(pts))

vertices(p::ConvexPolygon) = p.contour
nvertices(p::ConvexPolygon) = length(p.contour)
rings(p::ConvexPolygon) = [p.contour]
                                     
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



#function Meshes.normal(p::ConvexPolygon{2})
#    n = normal_(p)
#    return n / abs(n)
#end


area(p::ConvexPolygon) = norm(normalarea(p))

normal(p::ConvexPolygon) = normalize(normalarea(p))

"""
`pnpoly(poly, p)`

Determine if a point is inside a polygon. If a point is on the boundary of the polygon
this function might return `true` or `false` depending on floating point errors and
such.
"""
function pnpoly(p::SVec{2,T}, poly::ConvexPolygon{2,T}; atol=atolf(T)) where {T}

    v = vertices(poly)

    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = (v₂ - v₁)⋅(p - v₁)
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end

function pnpoly(p::SVec{3,T}, poly::ConvexPolygon{3,T}; atol=atolf(T)) where {T}

    
    v = vertices(poly)
    n = normal(poly)

    pl = Plane(v[1], n) # Plane

    # Check whether the point is in the plane of polygon
    if pnplane(p, pl, atol=atol)
        return false
    end

    
    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = (v₂ - v₁)⋅(p - v₁)
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end
                
Base.in(pt::SVec, p::ConvexPolygon) = pnpoly(pt, p)


function centroid(p::ConvexPolygon{Dim,T}) where {Dim,T}
    pts = vertices(p)

    p1 = pts[begin]  # Let's work with vectors
    npts = nvertices(p)
    A = zero(T)
    CA = zero(p1)
    for i in firstindex(pts)+1:lastindex(pts)-1
        u = pts[i]-pts[1]
        v = ustrip.(pts[i+1]-pts[1])
        Ai = norm(u × v) / 2  # Area of the triangle
        A += Ai
        # Compute the centroid of the triangle: mean of the coordinates
        Ci = ( p1 + pts[i] + pts[i+1]) ./ 3
        CA += Ai .* Ci
    end
    return SVec(CA./A) 
            
end



#=
function Meshes.simplexify(poly::ConvexPolygon)

    pts = vertices(poly)

    conn = [connect((1,i,i+1)) for i in firstindex(pts)+1:lastindex(pts)-1]
    return SimpleMesh(pts, conn)
  
    
end
=#
#discretize(poly::ConvexPolygon, method=nothing) = simplexify(poly)

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


