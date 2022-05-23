import GeometryBasics: AbstractPoint, Point, Point2, Point3

export PressureTap, xcoord, ycoord, zcoord, coordinates

struct PressureTap{T} <: AbstractPoint{3,T}
    "Coordinate of the pressure tap"
    point::Point{3,T}
    "Face to which the the pressure tap belongs"
    face::Int
    "Integer id of the pressure tap"
    id::Int
    "Is the pressure tap on the inside surface of the face?"
    isint::Bool
    "Does the pressure tap present problems?"
    hasprob::Bool
end

PressureTap(x,y,z; face=1, id=1, isint=false, hasprob=false) =
    PressureTap(Point3{Float64}(x,y,z),face,id,isint, hasprob)

PressureTap(p::Point{3,T}; face=1, id=1, isint=false, hasprob=false) where {T} =
    PressureTap(p,face,id,isint, hasprob)

#import Base.getindex
#getindex(p::PressureTap, i) = p.position[i]

xcoord(p::PressureTap) = p.point[1]
ycoord(p::PressureTap) = p.point[2]
zcoord(p::PressureTap) = p.point[3]
coordinates(p::PressureTap) = p.point



import Base.convert
convert(::Type{Point{3,T}}, p::PressureTap{T}) where {T} = p.point


