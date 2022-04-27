import GeometryBasics: Point, Point3, AbstractPoint

export PressureTap, xcoord, ycoord, zcoord

struct PressureTap{T} #<: AbstractPoint{Dim,T}
    "Coordinate of the pressure tap"
    position::Point3{T}
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
    PressureTap(Point3(x,y,z),face,id,isint, hasprob)

PressureTap(p::Point3{T}; face=1, id=1, isint=false, hasprob=false) where {T} =
    PressureTap(p,face,id,isint, hasprob)

#import Base.getindex
#getindex(p::PressureTap, i) = p.position[i]

xcoord(p::PressureTap) = p.position[1]
ycoord(p::PressureTap) = p.position[2]
zcoord(p::PressureTap) = p.position[3]


import Base.convert
convert(::Type{Point{3,T}}, p::PressureTap{T}) where {T} = p.position

