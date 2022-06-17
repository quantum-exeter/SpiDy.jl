## Noise type ##
abstract type Noise end

struct ClassicalNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::ClassicalNoise) = ω -> iszero(ω) ? zero(ω) : 2*n.T/ω

struct QuantumNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::QuantumNoise) = ω -> iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

struct BarkerNoise{TT<:Real} <: Noise
  T::TT
end

spectrum(n::BarkerNoise) = ω -> iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/T/2) - sign(ω))