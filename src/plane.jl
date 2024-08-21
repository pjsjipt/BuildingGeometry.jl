

struct Plane{Dim,T<:AbstractFloat}
    p::Point{Dim,T}
    n::Vec{Dim,T}
end

Plane(p,n) = Plane(p, n ./ norm(n))

normal(pl::Plane) = pl.n

basis(pl::Plane) = (householderbasis(pl.n)..., n)

Base.in(pt::Point, pl::Plane) = isapproxzero(normal(pl)⋅(p-pl.pt))

