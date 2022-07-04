## Coupling type ##
abstract type Coupling end

## Anisotropic coupling structure ##
struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling
    C::TT
    AnisoCoupling(C::TT) where {TT<:AbstractMatrix{T} where {T<:Real}} = size(C) == (3,3) ? new{TT}(C) : error("Coupling matrix must be 3x3")
end

Base.size(coupling::AnisoCoupling{TT}) where {TT<:AbstractMatrix{T} where {T<:Real}} = (3,3)
Base.getindex(coupling::AnisoCoupling{TT}, i::Int) where {TT<:AbstractMatrix{T} where {T<:Real}} = coupling.C[i]
Base.getindex(coupling::AnisoCoupling{TT}, I::Vararg{Int,N}) where {TT<:AbstractMatrix{T} where {T<:Real},N} = coupling.C[I...]

# x::Real
#            y::Real
#            OrderedPair(x,y) = x > y ? error("out of order") : new(x,y)

## Isotropic coupling structure ##
struct IsoCoupling{TT<:Real} <: Coupling
    C::TT
end

Base.size(coupling::IsoCoupling{TT}) where {TT<:Real} = (3,3)
Base.getindex(coupling::IsoCoupling{TT}, I::Vararg{Int,N}) where {TT<:Real,N} = I[1] == I[2] ? coupling.C : zero(TT)