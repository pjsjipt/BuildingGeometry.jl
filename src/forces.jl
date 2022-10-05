# Force calculation


"""
`addforcecontrib!(F, nodes, forces=(1,2,6), sgn=1, side=1, point=Point(0,0,0))`


The pressure on each face of a surface contributes to the force. This function
adds the contribution of a face characterized by nodes.

The contribution is assembled into a force matrix that when multiplied by
the pressure, results in the forces specified by components `forces`.

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

"""
`forcematrix(F, nodes, forces=(1,2,6), sgn=1, side=1, point=Point(0,0,0))`

Creates a force matrix from a [`BuildingSurface`](@ref) object. This function
initializes the matrix and calls [`addforcecontrib!`](@ref).

"""
function forcematrix(ncols::Integer, msh::BuildingSurface{T}, forces=(1,2,6);
                     sgn=1, side=1, point=Point{3,T}(0,0,0)) where {T}
    
    F = zeros(T, length(forces), ncols)
    return addforcecontrib!(F, msh, forces; sgn=sgn, side=side, point=point)
end


"""
`addforcecontrib!(F, msh, forces; interleaved=false, sgn=1, side=1,
                          point=Point{3,T}(0,0,0))`

Adds the force contribution from a list of meshes given by argument `msh`. The specific
forces calculated is given by argument `forces`.

### `interleaved`

Each section of the mesh (an element of vector `msh`) has different components of force, specified by arguemtn `forces`.

If `interleaved == false`, the forces are numbered in the order of force component.

On the other hand, if `interleaved == true`, the forces are numbered each section in equence.

If the components of force are (1,2,6) (Fx, Fy, Mz), when `interleaved == false`,
when the forces are calculated, they will be in the following order:

Fx₁
Fx₂
...
Fxₙ
Fy₁
Fy₂
...
Fyₙ
Mz₁
Mz₂
...
Mzₙ

In the list above, the subindices are the floor number (or index of `msh` vector).

When `interleaved == true`, the output forces will be ordered as:

Fx₁
Fy₁
Mz₁
Fx₂
Fy₂
Mz₂
...
Fxₙ
Fyₙ
Mzₙ

"""
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

"""                    
`forcematrix(F, msh, forces; interleaved=false, sgn=1, side=1,
                              point=Point{3,T}(0,0,0))`

Allocates memory and calls [`addforcecontrib!`](@ref) to compute the force matrix
for different sections (usually floors of a building).
"""                     
function forcematrix(ncols::Integer, msh::AbstractVector{<:BuildingSurface{T}},
                     forces=(1,2,6); interleaved=false, sgn=1, side=1,
                     point=Point{3,T}(0,0,0)) where {T}
    
    F = zeros(T, length(msh)*length(forces), ncols)
    return addforcecontrib!(F, msh, forces; interleaved=interleaved,
                            sgn=sgn, side=side, point=point)
end




                                            
