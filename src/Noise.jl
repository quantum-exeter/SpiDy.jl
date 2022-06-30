## Noise type ##
abstract type Noise end

## Classical noise structure ##
struct ClassicalNoise{TT<:Real} <: Noise
  T::TT
end

## Quantum noise structure ##
struct QuantumNoise{TT<:Real} <: Noise
  T::TT
end

## NoZeroQuantumNoise noise structure ##
struct NoZeroQuantumNoise{TT<:Real} <: Noise
  T::TT
end

"""
spectrum(n::ClassicalNoise)

Returns the classical noise at temperature n.T as a function of ω.
"""
spectrum(n::ClassicalNoise) = ω -> iszero(ω) ? zero(ω) : 2*n.T/ω

"""
spectrum(n::QuantumNoise)

Returns the quantum noise at temperature n.T as a function of ω.
"""
spectrum(n::QuantumNoise) = ω -> iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

"""
spectrum(n::NoZeroQuantumNoise)

Returns the quantum noise at temperature n.T as a function of ω. It differs from spectrum(n::QuantumNoise)
in the fact that the zero-point fluctuations are removed.
"""
spectrum(n::NoZeroQuantumNoise) = ω -> iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/T/2) - sign(ω))