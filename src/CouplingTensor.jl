"""
    abstract type Coupling

Abstract type used to represent different forms of the environment coupling tensor.
"""
abstract type Coupling end

"""
    struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling

A subtype of `Coupling` used to represent generic anisotropic coupling between the spin and environment.

# Fields
- `C::TT`: The coupling matrix, which must be a 3x3 matrix (of type `Matrix{<:Real}`).

# Examples
For a system with environment coupling where the `y` coupling is twice as large as the `x`,
and the `z` three times as large, one can do:
```julia
julia> C = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  2.0  0.0
 0.0  0.0  3.0

julia> coupling = AnisoCoupling(C)
AnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])

One can define effective 2D couplings by setting all coefficients in one of the dimensions to zero:
```julia
julia> C = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
3x3 Matrix{Float64}
 1.0 0.0 0.0
 0.0 1.0 0.0
 0.0 0.0 0.0

julia> coupling_2d = AnisoCoupling(C)
AnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0])
```
Similarly, one can define effective 1D couplings by setting two of the components to zero.
"""
struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling
    C::TT
    AnisoCoupling(C::TT) where {TT<:AbstractMatrix{T} where {T<:Real}} = size(C) == (3,3) ? new{TT}(C) : error("Coupling matrix must be 3x3")
end
Base.size(coupling::AnisoCoupling{TT}) where {TT<:AbstractMatrix{T} where {T<:Real}} = (3,3)
Base.getindex(coupling::AnisoCoupling{TT}, i::Int) where {TT<:AbstractMatrix{T} where {T<:Real}} = coupling.C[i]
Base.getindex(coupling::AnisoCoupling{TT}, I::Vararg{Int,N}) where {TT<:AbstractMatrix{T} where {T<:Real},N} = coupling.C[I...]

"""
    struct IsoCoupling{TT<:Real} <: Coupling

A subtype of `Coupling` used to represent isotropic coupling between spin and environment.

# Fields
- `C::Real`: The coupling strength.

# Examples
```julia
julia> coupling = IsoCoupling(2.5)
IsoCoupling{Float64}(2.5)
```
"""
struct IsoCoupling{TT<:Real} <: Coupling
    C::TT
end
Base.size(coupling::IsoCoupling{TT}) where {TT<:Real} = (1,1)
Base.getindex(coupling::IsoCoupling{TT}, I::Vararg{Int,N}) where {TT<:Real,N} = I[1] == I[2] ? coupling.C : zero(TT)