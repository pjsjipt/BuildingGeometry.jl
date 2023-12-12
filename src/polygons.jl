



struct ConvexPolygon{Dim,T,R<:Ring{Dim,T}} <: Polygon{Dim,T}
    contour::R
    function ConvexPolygon{Dim,T,R}(contour) where {Dim,T,R}
        @assert isclosed(contour)
        new(contour)
    end
end
"""
`ConvexPolygon(contour)`
`ConvexPolygon(x,y)`
`ConvexPolygon(x,y,z)`
`ConvexPolygon(pts)`

Create a convex polygon. 


"""
ConvexPolygon(contour::R) where {Dim,T,V,R<:Ring{Dim,T,V}} =
    ConvexPolygon{Dim,T,Ring{Dim,T,V}}(contour)


ConvexPolygon(x::AbstractVector, y::AbstractVector) =
    ConvexPolygon([Point2(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    ConvexPolygon([Point3(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

function ConvexPolygon(pts::AbstractVector{P}) where {P<:Point}
    ConvexPolygon(Ring(pts))
end

ConvexPolygon(pts::AbstractVector{TP}) where {TP<:Tuple} = ConvexPolygon(Point.(pts))
ConvexPolygon(pts::Vararg{P}) where {P<:Point} = ConvexPolygon(collect(pts))

Meshes.vertices(p::ConvexPolygon) = vertices(p.contour)
Meshes.nvertices(p::ConvexPolygon) = nvertices(p.contour)

                                     
"""
`normal_(p)`

Compute the area vector of a convex polygon. By area vector, it is meant the outer
 vector whose length is the area of the convex polygon.
"""
function normal_(p::ConvexPolygon{2,T}) where {T}
    pts = vertices(p)
    A = zero(T)
    for i in firstindex(pts)+1:lastindex(pts)
        A += pts[i-1].coords[1]*pts[i].coords[2] - pts[i].coords[1]*pts[i-1].coords[2]
    end
    return A/2
end

function Meshes.normal(p::ConvexPolygon{2})
    n = normal_(p)
    return n / abs(n)
end


Meshes.area(p::ConvexPolygon{2}) = abs(normal_(p))

"""
`pnpoly(poly, p)`

Determine if a point is inside a polygon. If a point is on the boundary of the polygon
this function might return `true` or `false` depending on floating point errors and
such.
"""
function pnpoly(poly::ConvexPolygon{2}, p::Point{2})

    v = vertices(poly)

    v₁ = v[begin]

    for i in eachindex(v)
        v₂ = v[i+1]
        s = (v₂ - v₁) ⋅ (p - v₁)
        if s < 0
            return false
        end
        v₁ = v₂
    end
    return true
    
end

Base.in(pt::Point, p::ConvexPolygon) = pnpoly(p, pt)


    
                                   
function normal_(p::ConvexPolygon{3,T}) where {T}
    v = vertices(p)
    nv = length(v)
    x = y = z = zero(T)
    u₁ = v[begin+1] - v[begin]
    for i in firstindex(v)+1:lastindex(v)-1
        u₂ = v[i+1]-v[i]
        prd = u₁ × u₂
        x += prd[1]
        y += prd[2]
        z += prd[3]
        u₁ = u₂
    end

    return Vec{3,T}(x/2, y/2, z/2)
    
end

function normal(p::ConvexPolygon{3})
    n = normal_(p)
    return n ./ norm(n)
end

    
Meshes.area(p::ConvexPolygon{3}) = hypot(normal_(p)...)

Meshes.measure(p::ConvexPolygon) = area(p)

"""
`centroid(p)`
Computes the centroid of a [`ConvexPolygon`](@ref).
"""                      
function Meshes.centroid(p::ConvexPolygon{2,T}) where {T}

    A = normal_(p)

    v = vertices(p)
    Cx = Cy = zero(T)
    
    for i in firstindex(v):lastindex(v)-1
        tmp = v[i].coords[1]*v[i+1].coords[2] - v[i+1].coords[1]*v[i].coords[2]
        Cx += (v[i].coords[1] + v[i+1].coords[1]) * tmp
        Cy += (v[i].coords[2] + v[i+1].coords[2]) * tmp
    end

    return Point{2,T}(Cx/6A, Cy/6A)
end

function centroid(p::ConvexPolygon{3,T}) where {T}

    A = zero(T)
    C = Vec{3,T}(0,0,0)
    v = vertices(p)
    
    v₀ = v[begin]
    for i = (firstindex(v)+1):(lastindex(v)-1)
        centr = (v₀.coords + v[i].coords + v[i+1].coords) ./ 3
        Aᵢ = norm( (v[i]-v₀) × (v[i+1]-v₀)) / 2
        A += Aᵢ
        C += Aᵢ * centr 
    end
    return Point{3,T}(C[1]/A, C[2]/A, C[3]/A)
    
end


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


