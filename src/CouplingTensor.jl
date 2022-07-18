"""
```Julia
Coupling
```

Definition of the abstract type `Coupling`.
"""
abstract type Coupling end

"""
```Julia
AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}}
```

Returns a `AnisoCoupling` structure of type `Coupling` built by passing it a 3x3 `Matrix{Real}` which defines
the n-dimentional coupling between the spin and the stochastic fields. The matrix can therefore
define a 1D, 2D as well as a 3D coupling. 
"""
struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling
    C::TT
    AnisoCoupling(C::TT) where {TT<:AbstractMatrix{T} where {T<:Real}} = size(C) == (3,3) ? new{TT}(C) : error("Coupling matrix must be 3x3")
end
Base.size(coupling::AnisoCoupling{TT}) where {TT<:AbstractMatrix{T} where {T<:Real}} = (3,3)
Base.getindex(coupling::AnisoCoupling{TT}, i::Int) where {TT<:AbstractMatrix{T} where {T<:Real}} = coupling.C[i]
Base.getindex(coupling::AnisoCoupling{TT}, I::Vararg{Int,N}) where {TT<:AbstractMatrix{T} where {T<:Real},N} = coupling.C[I...]

"""
```Julia
IsoCoupling{TT<:Real}
```

Returns a IsoCoupling of type `Coupling` built by passing it a single `Real` value which defines
the n-dimentional isotropic coupling between the spin and the stochastic fields.
"""
struct IsoCoupling{TT<:Real} <: Coupling
    C::TT
end
Base.size(coupling::IsoCoupling{TT}) where {TT<:Real} = (1,1)
Base.getindex(coupling::IsoCoupling{TT}, I::Vararg{Int,N}) where {TT<:Real,N} = I[1] == I[2] ? coupling.C : zero(TT)