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
function addforcecontrib!(F::AbstractMatrix{T}, msh::BuildingSurface{T}, forces=(1,2,6);
                          sgn=1, side=1, point=Point{3,T}(0,0,0)) where {T}
    nf = length(forces)
    (nf != size(F,1)) && error("`forces` and `F` have incompatible lengths")
    
    for n in msh.nodes
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
            F[kj,i] -= sgn*f⃗[k]
        end
    end
    return F
end

function forcematrix(ncols::Integer, msh::BuildingSurface{T}, forces=(1,2,6);
                     sgn=1, side=1, point=Point{3,T}(0,0,0)) where {T}
    
    F = zeros(T, length(forces), ncols)
    return addforcecontrib!(F, msh, forces; sgn=sgn, side=side, point=point)
end

                     
function addforcecontrib!(F::AbstractMatrix{T},
                          msh::AbstractVector{<:BuildingSurface{T}},
                          forces=(1,2,6); interleaved=false, sgn=1, side=1,
                          point=Point{3,T}(0,0,0)) where {T}
    
    nnodes = length(msh)

    nf = length(forces)
    (nnodes*nf != size(F,1)) && error("size(F,1) should equal to the number of meshes X number of forces ($(nnodes*nf))")

    
    if isa(point, Point)
        p = repeat([point], nnodes)
    elseif isa(point, AbstractVector{<:Point})
        p = point
    else
        error("Invalid point: $point of invalid type!")
    end
    npts = length(p)
    (npts != nnodes) && error("Number of points should the same as the number of nodes")

    
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
    
    for (i,m) in enumerate(msh)
        addforcecontrib!(view(F, rows[:,i], :), m, forces;
                         sgn=sgn, side=side, point=p[i])
    end
    return F
end

                    
                     
function forcematrix(ncols::Integer, msh::AbstractVector{<:BuildingSurface{T}},
                     forces=(1,2,6); interleaved=false, sgn=1, side=1,
                     point=Point{3,T}(0,0,0)) where {T}
    
    F = zeros(T, length(msh)*length(forces), ncols)
    return addforcecontrib!(F, msh, forces; interleaved=interleaved,
                            sgn=sgn, side=side, point=point)
end




                                            
