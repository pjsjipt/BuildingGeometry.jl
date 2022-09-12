# Force calculation


function addforcecontrib!(F, nodes::AbstractVector{<:NodeInfo{Dim,T}}, forces=(1,2,6);
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
end

                     
                     
                                            
