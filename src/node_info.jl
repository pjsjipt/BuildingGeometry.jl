
"""
`NodeInfo(A,iex,iin)`

Represents node information. In this context, node information
is information on outer and inner pressure tap.

In the simplest case, the external node is the id of the external pressure
tap and inner node the id of internal pressure tap. Often the internal
pressure is not measured and some other mean of calculation is used.

The pressure information can refer to an element of surface or a node.
So to calculate forces and moments, the contribution of this node is given
by the area times the normal
"""
struct NodeInfo{Dim,T,TSide}
    "Area times outward normal"
    A::Vec{Dim,T}
    "Coordinates of the node"
    point::Point{Dim,T}
    "Information on node sides"
    side::TSide
end

nodearea(n::NodeInfo) = n.A
nodepoint(n::NodeInfo) = n.point
extnode(n::NodeInfo) = n.side[1]
intnode(n::NodeInfo) = n.side[2]
nodeside(n::NodeInfo,i) = n.side[i]
nodeside(n::NodeInfo) = n.side


    
