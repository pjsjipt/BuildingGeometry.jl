import StaticArrays: SA, SMatrix
import Unitful
import Unitful: uparse, ustrip
"""
`reescalemesh(tri; scale=1/200, origin=SVec(0.,0.,0.), unit="mm")`

Creates a new mesh that is simply a scaled version of mesh `tri`.
 The function assumes that the input mesh is in model scale will generate
a mesh in prototype scale. This means that if `scale` < 1, the new mesh
will be larger.

The measurement units of input mesh `tri` is given by keyword argument `unit`.
The output unit is given by keyword argument `ounit`.

The origin can be changed using keyword argument `origin`.
This is a point *in the original coordinate system* that will serve as the
coordinate reference point in the output.

## Arguments
 * `geom` A geometry object to be scaled
 * `scale`: scale factor. Coordinates will be *divided* by this number
 * `origin`: Origin of the coordinate system of the input mesh in the output coordinate system
 * `unit`: A string containing the input mesh length unit system.
 * `ounit`: A string containing the output mesh length unit system.
"""
function reescalemesh(tri::Tri{Dim,T}; scale=T(1/200),
                      origin=SVec{Dim,T}(zero(T),0,0),
                      unit="mm", ounit="m") where {Dim,T}
    ufac = ustrip(uparse(ounit), one(T) * uparse(unit)) / scale
    return reescalemesh(tri, ufac, origin) 
end


function reescalemesh(tri::Tri{Dim,T}, ufac::T,
                      origin=SVec{Dim,T}(zero(T),0,0)) where {Dim,T}
    v = vertices(tri)
    Tri( ( (v[i] - origin) .* ufac for i in 1:nvertices(tri) )... )
end
function reescalemesh(p::SVec{Dim,T}; scale=T(1/200), origin=SVec{Dim,T}(zero(T),0,0),
                      unit="mm", ounit="m") where {Dim,T}
    ufac = ustrip(uparse(ounit), one(T) * uparse(unit)) / scale
    return (p-origin) .* ufac 
end

reescalemesh(p::SVec{Dim,T}, ufac::T,
             origin=SVec{Dim,T}(zero(T),0,0)) where {Dim,T} =
                 (p - origin) .* ufac

"""
`translatemesh(tri, u)`

Translate a geometric object by vector `u`
"""
translatemesh(tri::Tri{Dim,T}, u::SVec{Dim,T}) where {Dim,T} =
    Tri( (v + u for v in vertices(tri) )...)

translatemesh(p::SVec, u::SVec) = p + u

"""
`rotationmatrix(θ, ω=SVec(0.0, 0.0, 1.0))`

Calculate the rotation matrix around vector `ω` by an angle of `θ`.
"""
function rotationmatrix(θ::T, ω::SVec{3,T}) where {T}
    c = cos(θ)
    s = sin(θ)
    u,v,w = ω / norm(ω)
       
    SA[ (c + u*u*(1-c))  (u*v*(1-c)-w*s)  (u*w*(1-c)+v*s);
        (v*u*(1-c)+w*s)  (c+v*v*(1-c))    (v*w*(1-c)-u*s);
        (w*u*(1-c)-v*s)  (w*v*(1-c)+u*s)  (c + w*w*(1-c))]
end

function rotationmatrix(θ::T) where {T}
    c = cos(θ)
    s = sin(θ)
       
    SA[ c   -s;
        s  c]
end


"""
`rotatemesh(msh, θ, ω, p)`

Rotates a mesh (or geometry) by angle `θ` around an axis `ω` passing through point
`p0`.
"""
function rotatemesh(tri::Tri{3,T}, θ::T,
                    ω::SVec{3,T}, p=zero(ω)) where {T}
    R = rotationmatrix(θ, ω)
    return rotatemesh(tri, R, p)
end

    
rotatemesh(tri::Tri{Dim,T}, R::SMatrix{T},
           p=zero(vertices(tri)[begin])) where {Dim,T} =
               Tri( (R*(v - p) + p for v in vertices(geom) )...)


rotatemesh(pt::SVec{Dim,T}, R::SMatrix{T}, p=zero(pt)) where {Dim,T} =
    (R*(pt-p) + p)

rotatemesh(pt::SVec{3,T}, θ::T, ω::SVec{3,T}, p=zero(ω)) where {T} =
    rotatemesh(pt, rotationmatrix(θ, ω), p)

rotatemesh(tri::Tri{2,T}, θ::T, p=SVec2{T}(zero(T),0)) where {T} =
    roratemesh(tri, rotationmatrix(θ), p)


rotatemesh(pt::SVec2{T}, θ::T, p=zero(pt)) where {T} =
    rotatemesh(pt, rotationmatrix(θ), p)

    
