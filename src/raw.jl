
# reading raw files

"""
`readraw(fname)`

Reads a `.raw` 3d file into a list of surfaces discretized into triangles.
"""
function readraw(fname)
    TriFace = Tri{3,Float64}
    mshlst = Vector{TriFace}[]
    open(fname, "r") do io
        msh = TriFace[]
        newobject = false
        while !eof(io)
            line = lowercase(readline(io))
            if occursin(r"object[0-9]+", line)
                # New object
                if length(msh) > 0
                    push!(mshlst, msh)
                    msh = TriFace[]
                end
            else # It should be a triangle
                vals = parse.(Float64, split(line, (' ', '\t'), keepempty=false))
                if length(vals) != 9
                    error("Unable to read triangle from raw file: wrong length!")
                end
                push!(msh, Tri(SVec(vals[1], vals[2], vals[3]),
                               SVec(vals[4], vals[5], vals[6]),
                               SVec(vals[7], vals[8], vals[9])))
            end
        end
        push!(mshlst, msh)
    end
    return mshlst

end
