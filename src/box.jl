

struct Box{Dim,T} <: GeometryPrimitive{Dim,T}
    min::Point{Dim,T}
    max::Point{Dim,T}

end

Base.maximum(b::Box) = b.max
Base.minimum(b::Box) = b.min
Base.extrema(b::Box) = (b.min, b.max)

function Base.in(p::Point{Dim}, box::Box{Dim}) where {Dim}
    all(box.min[i] <= p[i] <= box.max[i] for i in eachindex(p))
end

GeometryBasics.area(b::Box{2}) = (b.max[1] - b.min[1]) * (b.max[2] - b.min[2])
GeometryBasics.volume(b::Box{3}) = (b.max[1] - b.min[1]) * (b.max[2] - b.min[2]) *
    (b.max[3] - b.min[3])

centroid(b::Box{Dim}) where {Dim} = Point( ((b.min[i]+b.max[i])/2 for i in 1:Dim)... )


