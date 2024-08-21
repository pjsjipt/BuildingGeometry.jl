# Functions takes from Meshes.jl
# Meshes.jl/src/utils/misc.jl


"""
    householderbasis(n)

Returns a pair of orthonormal tangent vectors `u` and `v` from a normal `n`,
such that `u`, `v`, and `n` form a right-hand orthogonal system.

## References

* D.S. Lopes et al. 2013. ["Tangent vectors to a 3-D surface normal: A geometric tool
  to find orthogonal vectors based on the Householder transformation"]
  (https://doi.org/10.1016/j.cad.2012.11.003)
"""
function householderbasis(n::Vec{3,T}) where {T<:AbstractFloat}
  n̂ = norm(n)
  i = argmax(n .+ n̂)
  n̂ᵢ = Vec(ntuple(j -> j == i ? n̂ : zero(T), 3))
  h = n + n̂ᵢ
  H = (I - 2h * transpose(h) / (transpose(h) * h)) 
  u, v = [H[:, j] for j in 1:3 if j != i]
  i == 2 && ((u, v) = (v, u))
  Vec(u), Vec(v)
end

