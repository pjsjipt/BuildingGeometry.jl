

struct Plane{Dim,T<:AbstractFloat}
    p::Point{Dim,T}
    n::Vec{Dim,T}
    Plane(p::Point{Dim,T}, n::Vec{Dim,T}) where {Dim,T} = new{Dim,T}(p, n ./ norm(n))
end


normal(pl::Plane) = pl.n

basis(pl::Plane) = (householderbasis(pl.n)..., n)

Base.in(pt::Point, pl::Plane; atol) = isapproxzero(normal(pl)⋅(p-pl.p))

