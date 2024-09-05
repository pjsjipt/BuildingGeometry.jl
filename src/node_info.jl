
"""
`NodeInfo(A,iex,iin)`

Represents node information. In this context, node information
is information on outer and inner pressure tap.

In the simplest case, the external node is the id of the external pressure
tap and inner node the id of internal pressure tap. Often the internal
pressure is not measured and some other mean of calculation is used.

The pressure information can refer to an element of surface or a node.
So to calculate forces and moments, the contribution of this node is given
by the area times the normal.

The `tag` field stores a tag for each side of the node. It is a tuple of
two integers that can be used to access different sections of the total mesh.

"""
struct NodeInfo{T,TSide}
    "Area times outward normal"
    A::SVector{3,T}
    "Coordinates of the node"
    point::SVector{3,T}
    "Information on node sides"
    side::TSide
    "Tag identifying each side"
    tag::Tuple{Int,Int}
end

nodeinfo(A,point,side,tag=(0,0)) = NodeInfo(ustrip.(SVector(A)),ustrip.(SVector(point)),side,(0,0))

function nodeinfo(A::Vec, point::Point, side, tag=(0,0))
    Ax = ustrip.(uconvert.(u"m", A) .* u"m")
    Px = ustrip.(uconvert.(u"m", to(point)))
    NodeInfo(Ax, Px, side, tag)
end


nodearea(n::NodeInfo) = n.A
nodepoint(n::NodeInfo) = n.point
extnode(n::NodeInfo) = n.side[1]
intnode(n::NodeInfo) = n.side[2]
nodeside(n::NodeInfo,i) = n.side[i]
nodeside(n::NodeInfo) = n.side
nodetag(n::NodeInfo) = n.tag
nodetag(n::NodeInfo,i) = n.tag[i]


    
