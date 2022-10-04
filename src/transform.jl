import Meshes: Geometry
import StaticArrays: SA, SMatrix
import Unitful
import Unitful: uparse, ustrip
"""
`reescalemesh(tri; scale=1/200, origin=Point(0,0,0), unit="mm")`

Creates a new mesh that is simply a scaled version of mesh `tri`.
 The function assumes that the input mesh is in model scale will generate
a mesh in prototype scale. This means that if `scale` < 1, the new mesh
will be larger.

The measurement units of input mesh `tri` is given by keyword argument `unit`.
The output unit is given by keyword argument `ounit`.

The origin can be changed using keyword argument `origin`.
This is a point *in the output coordinate system* that corresponds to the origin
of the input mesh coordinate system.

## Arguments
 * `geom` A geometry object to be scaled
 * `scale`: scale factor. Coordinates will be *divided* by this number
 * `origin`: Origin of the coordinate system of the input mesh in the output coordinate system
 * `unit`: A string containing the input mesh length unit system.
 * `ounit`: A string containing the output mesh length unit system.
"""
function reescalemesh(geom::G; scale=1/200,
                      origin=Point(0,0,0), unit="mm", ounit="m") where {G}
    ufac = ustrip(uparse(ounit), 1.0 * uparse(unit)) / scale
    return reescalemesh(geom, ufac, origin) 
end

reescalemesh(geom::G, ufac, origin=Point(0,0,0)) where {G} =
    G((vertices(geom) .- Point(0,0,0)) .* ufac .+ origin )

function reescalemesh(p::Point; scale=1/200, origin=Point(0,0,0),
                      unit="mm", ounit="m")
    ufac = ustrip(uparse(ounit), 1.0 * uparse(unit)) / scale
    return (p-Point(0,0,0,)) * ufac + origin
end

reescalemesh(p::Point, ufac, origin=Point(0,0,0)) = (p-Point(0,0,0))*ufac + origin

"""
`translatemesh(tri, u)`

Translate a geometric object by vector `u`
"""
translatemesh(geom::G, u::Vec) where {G<:Geometry} =
    G([v + u for v in vertices(geom)])

translatemesh(p::Point, u::Vec) = p + u

"""
`rotationmatrix(θ, ω=Vec(0.0, 0.0, 1.0))`

Calculate the rotation matrix around vector `ω` by an angle of `θ`.
"""
function rotationmatrix(θ, ω=Vec(0.0, 0.0, 1.0))
    c = cos(θ)
    s = sin(θ)
    u,v,w = ω / norm(ω)
       
    SA[ (c + u*u*(1-c))  (u*v*(1-c)-w*s)  (u*w*(1-c)+v*s);
        (v*u*(1-c)+w*s)  (c+v*v*(1-c))    (v*w*(1-c)-u*s);
        (w*u*(1-c)-v*s)  (w*v*(1-c)+u*s)  (c + w*w*(1-c))]
end

"""
`rotatemesh(msh, θ, ω, p0)`

Rotates a mesh (or geometry) by angle `θ` around an axis `ω` passing through point
`p0`.
"""
function rotatemesh(geom::G, θ, ω::Vec, p0=Point(0,0,0)) where {G<:Geometry}
    R = rotationmatrix(θ, ω)
    
    return rotatemesh(geom, R, p0)
end

rotatemesh(geom::G, R::SMatrix, p0=Point(0,0,0)) where {G<:Geometry} =
    G([R*(v - p0) + p0 for v in vertices(geom)])

rotatemesh(p::Point, R::SMatrix, p0=Point(0,0,0)) = (R*(p-p0) + p0)
rotatemesh(p::Point, θ::Real, ω::Vec, p0=Point(0,0,0)) =
    rotatemesh(p, rotationmatrix(θ, ω), p0)
