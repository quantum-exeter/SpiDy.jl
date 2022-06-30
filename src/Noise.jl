## Noise type ##
abstract type Noise end

## Classical noise ##
struct ClassicalNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::ClassicalNoise) = ω -> iszero(ω) ? zero(ω) : 2*n.T/ω

## Quantum noise ##
struct QuantumNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::QuantumNoise) = ω -> iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

## NoZeroQuantumNoise noise ##
struct NoZeroQuantumNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::NoZeroQuantumNoise) = ω -> iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/T/2) - sign(ω))