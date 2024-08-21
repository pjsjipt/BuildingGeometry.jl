
import CircularArrays: CircularVector



struct ConvexPolygon{Dim,T<:AbstractFloat} <: AbstractPolygon{Dim,T}
    contour::CircularVector{Point{Dim,T}}
end
"""
`ConvexPolygon(contour)`
`ConvexPolygon(x,y)`
`ConvexPolygon(x,y,z)`
`ConvexPolygon(pts)`

Create a convex polygon. 


"""
function ConvexPolygon(contour::AbstractVector{Point{Dim,T}}) where {Dim,T}
    ConvexPolygon(CircularVector(contour))
end


function ConvexPolygon(x::AbstractVector{U}, y::AbstractVector{V}) where{U,V}
    T = promote_type(U,V)
    if T <: Integer
        T = Float64
    end
    ConvexPolygon([Point{2,T}(T(xx), T(yy)) for (xx,yy) in zip(x,y)])
end

function ConvexPolygon(x::AbstractVector{U}, y::AbstractVector{V},
                       z::AbstractVector{W}) where{U,V,W}
    T = promote_type(U,V,W)
    if T <: Integer
        T = Float64
    end
    ConvexPolygon([Point{3,T}(T(xx), T(yy), T(zz)) for (xx,yy,zz) in zip(x,y,z)])
end

function ConvexPolygon(pts::AbstractVector{NTuple{N,T}}) where {N,T}
    T1 = (T <: AbstractFloat) ? T : Float64
    ConvexPolygon([Point( (T1.(p))... ) for p in pts])
end

                   

ConvexPolygon(pts...) = ConvexPolygon(collect(pts))

GeometryBasics.coordinates(p::ConvexPolygon) = p.contour
Base.length(p::ConvexPolygon) = length(p.contour)



"""
`normalarea(p)`

Compute the area vector of a convex polygon. By area vector, it is meant the outer
 vector whose length is the area of the convex polygon.
"""
normalarea(pts::AbstractVector{Point{Dim,T}}) where {Dim,T} =
    sum( (pts[i-1]-pts[begin]) × (pts[i]-pts[begin])
         for i in firstindex(pts)+1:lastindex(pts) ) / 2

normalarea(p::ConvexPolygon) = normalarea(coordinates(p))

GeometryBasics.area(p::ConvexPolygon) = norm(normalarea(p))

function normal(p::ConvexPolygon)
    na = normalarea(p)

    return na ./ norm(na)
end

    



"""
`pnpoly(poly, p)`

Determine if a point is inside a polygon. If a point is on the boundary of the polygon
this function might return `true` or `false` depending on floating point errors and
such.
"""
function pnpoly(poly::ConvexPolygon{2}, p::Point{2})

    v = coordinates(poly)

    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = (v₂ - v₁) ⋅ (p - v₁) 
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true

    
end


function pnpoly(poly::ConvexPolygon{3}, p::Point{3})

    v = coordinates(poly)
    n = normal(poly)
    
    # Check
    pl = Plane(v[1], n) # Plane

    # Check whether the point is in the plane of polygon
    if p ∉ pl
        return false
    end

    
    v₁ = v[begin]

    for v₂ in v[begin+1:end+1]
        s = (v₂ - v₁) ⋅ (p - v₁)
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end


Base.in(pt::Point, p::ConvexPolygon) = pnpoly(p, pt)


"""
`centroid(p)`
Computes the centroid of a [`ConvexPolygon`](@ref).
"""                      
function centroid(p::ConvexPolygon{Dim,T}) where {Dim,T}
    pts = coordinates(p)

    p1 = pts[begin]  # Let's work with vectors
    N = Dim
    A = zero(T)
    CA = zero(p1)
    for i in firstindex(pts)+1:lastindex(pts)-1
        u = pts[i]-pts[1]
        v = pts[i+1]-pts[1]
        
        Ai = norm(u × v) / 2  # Area of the triangle
        A += Ai
        # Compute the centroid of the triangle: mean of the coordinates
        Ci = ( p1 + pts[i] + pts[i+1] ) ./ 3
        CA += Ai .* Ci
    end
    return  CA./A
    
end





"""
`poly2mesh(poly)`

Creates a triangle mesh from a convex polygon. 
"""
function poly2mesh(poly::ConvexPolygon)
    pts = coordinates(poly)
    
   
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



