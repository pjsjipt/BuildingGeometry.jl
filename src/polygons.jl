

import GeometryBasics: AbstractPoint, Point, Point2, AbstractPolygon, area, Point3, Polygon

import GeneralPolygonClipper as gpc
  


export SimplePolygon, area, centroid, normal, coordinates, simple2poly


struct SimplePolygon{Dim,T,P<:AbstractPoint{Dim,T},L<:AbstractVector{P}} <: AbstractPolygon{Dim,T}
    contour::L
end

SimplePolygon(x::AbstractVector, y::AbstractVector) =
    SimplePolygon([Point2(xx, yy) for (xx,yy) in zip(x,y)])

SimplePolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    SimplePolygon([Point3(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

#SimplePolygon(pts::E) where {E<:AbstractVector{P}} where {P<:AbstractPoint{Dim,T}} where {Dim,T} = SimplePolygon{Dim,T,P,E}(pts)

import GeometryBasics.coordinates
coordinates(p::SimplePolygon) = p.contour

area(p::SimplePolygon{2,T}) where {T} = area(coordinates(p))

normal(p::SimplePolygon{2,T}) where {T} = area(coordinages(p))

crossprod(u::Point3, v::Point3) = (u[2]*v[3] - u[3]*v[2],
                                   u[3]*v[1] - u[1]*v[3],
                                   u[1]*v[2] - u[2]*v[1])
                                   
function normal(p::SimplePolygon{3,T}) where {T}
    v = coordinates(p)
    nv = length(v)
    x,y,z = crossprod(v[end], v[begin])

    for i in 2:nv
        prd = crossprod(v[i-1], v[i])
        x += prd[1]
        y += prd[2]
        z += prd[3]
    end

    return Point3(x/2, y/2, z/2)
    
end

area(p::SimplePolygon{3}) = hypot(normal(p)...)

"""
`centroid(p)`

Computes the centroid of a [`SimplePolygon`](@ref).
"""                      
function centroid(p::SimplePolygon{2,T}) where {T}

    A = area(p)

    v = coordinates(p)
    nv = length(v)

    # xᵢyᵢ₊₁ - xᵢ₊₁yᵢ  First index = last index
    tmp = v[end][1]*v[begin][2] - v[begin][1]*v[end][2]
    
    Cx = (v[end][1] + v[begin][1]) * tmp
    Cy = (v[end][2] + v[begin][2]) * tmp

    for i in firstindex(v):lastindex(v)-1
        tmp = v[i][1]*v[i+1][2] - v[i+1][1]*v[i][2]
        Cx += (v[i][1] + v[i+1][1]) * tmp
        Cy += (v[i][2] + v[i+1][2]) * tmp
    end

    return Point2{T}(Cx/6A, Cy/6A)
end

function centroid(p::SimplePolygon{3,T}) where {T}
    v = coordinates(p)
    nv = length(v)
    nrm = normal(p)
    A = hypot(nrm...)
    n⃗ = nrm ./ A
    
    C = Point3{T}(0,0,0)
    v₀ = v[begin]
    for i = (firstindex(v)+1):(lastindex(v)-1)
        centr = (v₀ + v[i] + v[i+1]) ./ 3
        Aᵢ⃗ = crossprod(v[i]-v₀, v[i+1]-v₀) ./ 2
        s = sign(sum(Aᵢ⃗ .* n⃗))
        C += s * centr .* hypot(Aᵢ⃗...)
    end
    return C
    
end

import GeneralPolygonClipper: Vertex
import Base.union

"""
`gpc_operations(op::GPCOperation, p1, p2)`

Performs boolean operations on 2D polygons. First the polygons
are converted to [`GPCPolygon`](@ref). Then, performs the boolean
operation specified by parameter `op` (see enum [`GPCOperation`](@ref)).

Finally, convert the results to [`SimplePolygon`](@ref). There can be 
more than one resulting polygon depending on the exact nature of `p1`, `p2`
and `op`.

"""
function gpc_operations(op::gpc.GPCOperation, p1::SimplePolygon{2,Float64},
                        p2::SimplePolygon{2,Float64})

    gpc1 = gpc.GPCPolygon([false], [Vertex(p[1], p[2]) for p in coordinates(p1)])
    gpc2 = gpc.GPCPolygon([false], [Vertex(p[1], p[2]) for p in coordinates(p2)])

    gpcout = gpc.gpc_polygon_clip(op, gpc1, gpc2)

    # If both polygons were convex, the output would be a single polygon
    # In other cases there can be more than one output polygon.

    return SimplePolygon.([[Point2{Float64}(v.x, v.y) for v in g] for g in gpcout.contours])
    
end

import Base.union
union(p1::SimplePolygon{2,Float64}, p2::SimplePolygon{2,Float64}) =
    gpc_operations(gpc.GPC_UNION, p1, p2)
import Base.intersect
intersect(p1::SimplePolygon{2,Float64}, p2::SimplePolygon{2,Float64}) =
    gpc_operations(gpc.GPC_INT, p1, p2)
import Base.(-)
(-)(p1::SimplePolygon{2,Float64}, p2::SimplePolygon{2,Float64}) =
    gpc_operations(gpc.GPC_DIFF, p1, p2)

import Base.xor
xor(p1::SimplePolygon{2,Float64}, p2::SimplePolygon{2,Float64}) =
    gpc_operations(gpc.GPC_XOR, p1, p2)

simple2poly(p::SimplePolygon{2,Float64}) = Polygon(coordinates(p))
