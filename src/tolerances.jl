# Taken from Meshes.jl
# Meshes.jl/src/tolerances.jl
# Meshes.jl/src/utils/cmp.jl


const ATOL64 = 1.0e-10
const ATOL32 = 1.0f-5

atol(::Type{Float64}) = ATOL64
atol(::Type{Float32}) = ATOL32
atol(x) = atol(typeof(x))

isequalzero(x) = x == zero(x)
isequalone(x) = x == one(x)

isapproxequal(x, y; atol=atol(x), kwargs...) = isapprox(x, y; atol, kwargs...)
isapproxzero(x; atol=atol(x), kwargs...) = isapprox(x, zero(x); atol, kwargs...)
isapproxone(x; atol=atol(x), kwargs...) = isapprox(x, oneunit(x); atol, kwargs...)

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonpositive(x) = x ≤ zero(x)
isnonnegative(x) = x ≥ zero(x)
