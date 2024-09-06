
# Basic definitions


atolf(L::T, rtol::T) where {T}= L*rtol
atolf(::Type{Float64}) = 1e-10
atolf(::Type{Float32}) = 1e-5
atolf(L::T) where{T} = atolf(L, atolf(T))


const SVec{Dim,T<:AbstractFloat} = SVector{Dim,T}
const SVec2{T} = SVec{2,T}
const SVec3{T} = SVec{3,T}

struct Plane{Dim,T} <: AbstractBuildGeom
    p::SVec{Dim,T}
    n::SVec{Dim,T}
end


normal(p::Plane) = p.n

function pnplane(p::SVec{Dim,T}, plane::Plane{Dim,T}; atol=atolf(T)) where {Dim,T}
    abs( (p-plane.p)⋅normal(plane) ) < atol
end

Base.in(p::SVec{Dim,T}, plane::Plane{Dim,T}) where {Dim,T} = pnplane(p,plane)
  

struct Tri{Dim,T} <: AbstractBuildGeom
    v::SVector{3,SVec{Dim,T}}
end

Tri(v1,v2,v3) = Tri(@SVector [SVec(v1),SVec(v2),SVec(v3)])
Tri(v::AbstractVector{<:SVec}) = Tri(@SVector [v[1], v[2], v[3]])

nvertices(t::Tri) = 3
vertices(t::Tri) = t.v
normalarea(t::Tri) = (t.v[2]-t.v[1]) × (t.v[3]-t.v[1]) ./ 2
area(t::Tri) = norm(normalarea(t))
normal(t::Tri) = normalize(normalarea(t))
centroid(t::Tri) = sum(t.v) ./ 3

struct Box{Dim,T} <: AbstractBuildGeom
    min::SVec{Dim,T}
    max::SVec{Dim,T}
end


Base.extrema(b::Box) = (b.min,b.max)
Base.extrema(g::AbstractBuildGeom) = extrema(vertices(g))

function Base.extrema(pts::NTuple{2,SVec{Dim,T}}) where {Dim,T}
    pmin = ntuple(i->min(pts[1][i], pts[2][i]), Dim)
    pmax = ntuple(i->max(pts[1][i], pts[2][i]), Dim)
    SVec(pmin),SVec(pmax)
end

    
function Base.extrema(pts::AbstractVector{SVec{Dim,T}}) where {Dim,T}
    pmin = ntuple(i->minimum(p[i] for p in pts), Dim)
    pmax = ntuple(i->maximum(p[i] for p in pts), Dim)
    SVec(pmin),SVec(pmax)
end

function Base.extrema(geos::AbstractVector{<:AbstractBuildGeom})
    pmin,pmax = extrema(geos[begin])

    for i in firstindex(geos)+1:lastindex(geos)
        pmin1, pmax1 = extrema(geos[i])
        pmin, tmp = extrema((pmin,pmin1))
        pmax, tmp = extrema((pmax,pmax1))
    end
    SVec(pmin),SVec(pmax)
end

boundingbox(pts::AbstractVector{SVec{Dim,T}}) where {Dim,T} = Box(extrema(pts)...)
boundingbox(g::AbstractBuildGeom) = Box(extrema(g)...)

boundingbox(geos::AbstractVector{<:AbstractBuildGeom}) = Box(extrema(geos)...)




    
