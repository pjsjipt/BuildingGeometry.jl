

import GeometryBasics: AbstractPoint, Point, Point2, AbstractPolygon, area, Point3

export SimplePolygon, area, centroid, normal, coordinates

struct SimplePolygon{Dim,T,P<:AbstractPoint{Dim,T},L<:AbstractVector{P}} <: AbstractPolygon{Dim,T}
    contour::L
end


SimplePolygon(pts::E) where {E<:AbstractVector{P}} where {P<:AbstractPoint{Dim,T}} where {Dim,T} = SimplePolygon{Dim,T,P,E}(pts)

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
