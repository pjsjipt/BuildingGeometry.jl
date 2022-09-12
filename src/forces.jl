# Force calculation


"""
`addforcecontrib!(F, nodes, forces=(1,2,6), sgn=1, side=1, point=Point(0,0,0))`

Adds the force contribution of a face characterized by nodes.
Only the components in `forces` is used where the components are

 1. `Fx`
 2. `Fy`
 3. `Fz`
 4. `Mx`
 5. `My`
 6. `Mz`

The contribution of each triangle in the mesh should be along the outer normal.
To invert this, set `sgn=-1`. The side of the node is given by parameter `side`.
Nodes where the side index are 0 or negative, no direct measurements for this node
are available and its contribution should be added directly.

The moments are calculated with respect to point `point`.
"""
function addforcecontrib!(F::AbstractMatrix{T}, nodes::AbstractVector{<:NodeInfo{Dim,T}}, forces=(1,2,6);
                          sgn=1, side=1, point=Point{Dim,T}(0,0,0)) where {Dim,T}
    nf = length(forces)
    (nf != size(F,1)) && error("`forces` and `F` have incompatible lengths")
    
    for n in nodes
        i =  nodeside(n,side)
        if i <= 0
            # No pressure tap corresponding to it
            continue
        end
        A⃗ = nodearea(n)
        r⃗ = nodepoint(n) - point
        M⃗ = r⃗ × A⃗
        f⃗ = [A⃗; M⃗]
        for (kj,k) in enumerate(forces)
            F[kj,i] += sgn*f⃗[k]
        end
    end
    return F
end

function forcematrix(ncols::Integer, nodes::AbstractVector{<:NodeInfo{Dim,T}}, forces=(1,2,6);
                     sgn=1, side=1, point=Point{Dim,T}(0,0,0)) where {Dim,T}
    
    F = zeros(T, length(forces), ncols)
    return addforcecontrib!(F, nodes, forces; sgn=sgn, side=side, point=point)
end

                     
function addforcecontrib!(F::AbstractMatrix{T}, nodes::AbstractVector{MESH}, forces=(1,2,6); interleaved=false,
                          sgn=1, side=1, point::AbstractVector{Point{Dim,T}}) where {Dim,T,MESH<:AbstractVector{<:NodeInfo{Dim,T}}}
    npts = length(point)
    nnodes = length(nodes)

    (npts != nnodes) && error("Number of points should the same as the number of nodes")
    nf = length(forces)
    (nnodes*nf != size(F,1)) && error("size(F,1) should equal to the number of meshes X number of forces ($(nnodes*nf))")

    rows = zeros(Int,nf,nnodes)
    if interleaved
        cnt = 1
        for k in 1:nnodes
            for i in 1:nf
                rows[i,k] = cnt
                cnt += 1
            end
        end
    else
        cnt = 1 
        for i in 1:nf
            for k in 1:nnodes
                rows[i,k] = cnt
                cnt += 1
            end
        end
    end
    
    
    for (i,n) in enumerate(nodes)
        addforcecontrib!(view(F, rows[:,i], :), nodes[i], forces; sgn=sgn, side=side, point=point[i])
    end
    return F
end

                     
function addforcecontrib!(F::AbstractMatrix{T}, nodes::AbstractVector{MESH}, forces=(1,2,6); interleaved=false,
                          sgn=1, side=1, point=Point{Dim,T}(0,0,0)) where {Dim,T,MESH<:AbstractVector{<:NodeInfo{Dim,T}}}
    points = repeat(point, length(nodes))
    return addforcecontrib!(F, nodes, forces; interleaved=interleaved, sgn=sng, side=side, point=points)
end


                    
                     
function forcematrix(ncols::Integer, nodes::AbstractVector{MESH}, forces=(1,2,6); interleaved=false,
                     sgn=1, side=1, point::AbstractVector{Point{Dim,T}}) where {Dim,T,MESH<:AbstractVector{<:NodeInfo{Dim,T}}}
    F = zeros(T, length(nodes)*length(forces), ncols)
    return addforcecontrib!(F, nodes, forces; interleaved=interleaved, sgn=sgn, side=side, point=point)
end

function forcematrix(ncols::Integer, nodes::AbstractVector{MESH}, forces=(1,2,6); interleaved=false,
                     sgn=1, side=1, point=Point{Dim,T}(0,0,0)) where {Dim,T,MESH<:AbstractVector{<:NodeInfo{Dim,T}}}
    points = repeat(point, length(nodes))
    return forcematrix(ncols, nodes, forces; interleaved=interleaved, sgn=sgn, side=side, point=points)
end




                                            
