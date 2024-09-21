
# Data structures to project pressure measurement into surface data

export CpSurface, assemble!, assemble


struct CpSurface{T}
    "Assembly matrix for each mesh node"
    A::Matrix{T}
    "Degrees of freedom for each node"
    idx::Matrix{Int}
end

"""
`CpSurface(msh, side)`

Build a `CpSurface` object.

The `CpSurface` object projects the measurement array on each node of the mesh.

 * `msh::NodeInfo{T,Tuple{Int,Int}}`  Discretization of the building surface.
 * `side` side of the surface

The pressure on each node of the surface mesh corresponds to weighed average
of the pressures measured. This structure stores, for each node in the mesh,
the weights and channels.

For node `i`,
``pmsh[i] = sum(A[i,k]*p[idx[k]] for k in size(A,2))``

where `p` is the measured pressure and `pmsh` is the pressure on each node.

For now, only the simples cast is implemented: Each node is a triangle and
 corresponds to a single pressure, that is  `A[i,k] == 1` and `size(A,2) == 2`.

If no pressure tap is located on the triangle, just set A to zero and put a valid
pressure tap on idx.
"""
function CpSurface(msh::AbstractVector{NodeInfo{T,Tuple{Int,Int}}}, side) where {T}
    ntri = length(msh)
    idx = nodeside.(msh, side)
    A = ones(T,ntri,1)
    # Let's check the nodes where `idx < 0`  (nore pressure tap). We will set
    # idx to 1 (We know that this will be available...) and A = 0
    notap = idx .< 0
    A[notap,1] .= zero(T)
    idx[notap] .= 1
    return CpSurface(A,hcat(idx))
end

function assemble!(cpout::AbstractVector{T},
                   sdata::CpSurface{T}, cp::AbstractVector{T}) where {T}

    nout = length(cpout)
    @assert nout == size(sdata.A,1) "`cpout` should have one column per mesh node"

    for i in 1:nout
        cpout[i] = sdata.A[i] * cp[sdata.idx[i]]
    end

    return cpout
end


assemble(sdata::CpSurface{T}, cp::AbstractVector{T}) where {T} =
    assemble!(zeros(T, length(sdata.idx)), sdata, cp)


    
function assemble!(cpout::AbstractArray{T},
                   sdata::CpSurface{T}, cp::AbstractArray{T}) where {T}
    
    @assert size(cpout,1) == size(cp,1) "`cpout` should have the same number of rows as `cp`"
    @assert size(cpout,2) == size(sdata.A,1) "`cpout` should have one column per mesh node"
    # For now only the simplest case is implemented
    A = sdata.A
    idx = sdata.idx

    nt = size(cp,1)
    nout = size(cpout, 2)
    
    for (cpo, cpc) in zip(eachrow(cpout), eachrow(cp))
        for i in 1:nout
            cpo[i] = A[i] * cpc[idx[i]]
        end
    end
    return cpout
end

"""
`assemble(sdata, cp)`

Calculates the pressure from pressure measurement on each node of the mesh.
"""
assemble(sdata::CpSurface{T}, cp::AbstractArray{T}) where {T} =
    assemble!(zeros(T, size(cp,1), size(sdata.idx,1)), sdata, cp)

(sdata::CpSurface)(cp) = assemble(sdata, cp)
(sdata::CpSurface)(cpout, cp) = assemble!(cpout, sdata, cp)

    
