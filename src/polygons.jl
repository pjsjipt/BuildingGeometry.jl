

import GeometryBasics: AbstractPoint, Point, Point2, AbstractPolygon, area, Point3, Polygon

import GeneralPolygonClipper as gpc
  


export ConvexPolygon, area, centroid, normal, coordinates, simple2poly, pnpoly


struct ConvexPolygon{Dim,T,P<:AbstractPoint{Dim,T},L<:AbstractVector{P}} <: AbstractPolygon{Dim,T}
    contour::L
end

ConvexPolygon(x::AbstractVector, y::AbstractVector) =
    ConvexPolygon([Point2(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    ConvexPolygon([Point3(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

#ConvexPolygon(pts::E) where {E<:AbstractVector{P}} where {P<:AbstractPoint{Dim,T}} where {Dim,T} = ConvexPolygon{Dim,T,P,E}(pts)

import GeometryBasics.coordinates
coordinates(p::ConvexPolygon) = p.contour

function normal(p::ConvexPolygon{2})
    pts = coordinates(p)
    
    A = pts[end][1] * pts[begin][2] - pts[begin][1] * pts[end][2]

    for i in firstindex(pts)+1:lastindex(pts)
        A += pts[i-1][1]*pts[i][2] - pts[i][1]*pts[i-1][2]
    end
    return A/2
end


area(p::ConvexPolygon{2}) = abs(normal(p))
    
function pnpoly(poly::ConvexPolygon{2}, p::Point{2})

    test = false
    v = coordinates(poly)
    j = lastindex(v)
    for i in eachindex(v)
        if (  (v[i][2] > p[2]) != (v[j][2]>p[2]) ) &&
            ( p[1] ≤ (v[j][1]-v[i][1])*(p[2]-v[i][2]) / (v[j][2]-v[i][2]) + v[i][1])
            test = !test
        end
    end
    return test
    
end

Base.in(pt, p::ConvexPolygon) = pnpoly(p, pt)


    
crossprod(u::Point3, v::Point3) = (u[2]*v[3] - u[3]*v[2],
                                   u[3]*v[1] - u[1]*v[3],
                                   u[1]*v[2] - u[2]*v[1])
                                   
function normal(p::ConvexPolygon{3,T}) where {T}
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

area(p::ConvexPolygon{3}) = hypot(normal(p)...)



"""
`centroid(p)`
Computes the centroid of a [`ConvexPolygon`](@ref).
"""                      
function centroid(p::ConvexPolygon{2,T}) where {T}

    A = normal(p)

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

function centroid(p::ConvexPolygon{3,T}) where {T}
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
    return C/A
    
end


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



import GeneralPolygonClipper: Vertex
import Base.union

"""
`gpc_operations(op::GPCOperation, p1, p2)`
Performs boolean operations on 2D polygons. First the polygons
are converted to [`GPCPolygon`](@ref). Then, performs the boolean
operation specified by parameter `op` (see enum [`GPCOperation`](@ref)).
Finally, convert the results to [`ConvexPolygon`](@ref). There can be 
more than one resulting polygon depending on the exact nature of `p1`, `p2`
and `op`.
"""
function gpc_operations(op::gpc.GPCOperation, p1::ConvexPolygon{2,Float64},
                        p2::ConvexPolygon{2,Float64})

    gpc1 = gpc.GPCPolygon([false], [Vertex(p[1], p[2]) for p in coordinates(p1)])
    gpc2 = gpc.GPCPolygon([false], [Vertex(p[1], p[2]) for p in coordinates(p2)])

    gpcout = gpc.gpc_polygon_clip(op, gpc1, gpc2)

    # If both polygons were convex, the output would be a single polygon
    # In other cases there can be more than one output polygon.

    return ConvexPolygon.([[Point2{Float64}(v.x, v.y) for v in g] for g in gpcout.contours])
    
end

import Base.union
union(p1::ConvexPolygon{2,Float64}, p2::ConvexPolygon{2,Float64}) =
    gpc_operations(gpc.GPC_UNION, p1, p2)
import Base.intersect
intersect(p1::ConvexPolygon{2,Float64}, p2::ConvexPolygon{2,Float64}) =
    gpc_operations(gpc.GPC_INT, p1, p2)
import Base.(-)
(-)(p1::ConvexPolygon{2,Float64}, p2::ConvexPolygon{2,Float64}) =
    gpc_operations(gpc.GPC_DIFF, p1, p2)

import Base.xor
xor(p1::ConvexPolygon{2,Float64}, p2::ConvexPolygon{2,Float64}) =
    gpc_operations(gpc.GPC_XOR, p1, p2)

#simple2poly(p::SimplePolygon{2,Float64}) = Polygon(coordinates(p))
