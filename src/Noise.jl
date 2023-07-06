"""
    abstract type Noise

An abstract type used to represent different kinds of environment noise.

Any user-defined subtype of `Noise` should implement a method `noise(::NoiseSubType)`
which should return a function that represent the spectrum of the noise.
"""
abstract type Noise end

"""
    struct ClassicalNoise{TT<:Real} <: Noise

A subtype of `Noise` used to represent classical noise, with spectrum
```math
\\mathcal{S}(\\omega) = \\frac{k_\\mathrm{B}T}{\\hbar\\omega}
```

# Fields
- `T::TT`: The temperature of the environment associated with the noise.
"""
struct ClassicalNoise{TT<:Real} <: Noise
    T::TT
end

"""
    struct QuantumNoise{TT<:Real} <: Noise

A subtype of `Noise` used to represent quantum noise, with spectrum
```math
\\mathcal{S}(\\omega) = \\coth\\left(\\frac{\\hbar\\omega}{2k_\\mathrm{B}T}\\right)
```

# Fields
- `T::TT`: The temperature of the environment associated with the noise.
"""
struct QuantumNoise{TT<:Real} <: Noise
    T::TT
end

"""
    struct NoZeroQuantumNoise{TT<:Real} <: Noise

A subtype of `Noise` used to represent quantum noise with zero-point-fluctuations removed,
that is with spectrum
```math
\\mathcal{S}(\\omega) = \\coth\\left(\\frac{\\hbar\\omega}{2k_\\mathrm{B}T}\\right) - 1
```

# Fields
- `T::TT`: The temperature of the environment associated with the noise.
"""
struct NoZeroQuantumNoise{TT<:Real} <: Noise
    T::TT
end

# Returns the classical noise at temperature `n.T` as a function of `ω`. The conditional
# statement makes sure the function does not divide by zero in case of `ω==0`.
"""
    spectrum(n::ClassicalNoise)

Construct the spectrum of [`ClassicalNoise`](@ref) at a given frequency.

# Arguments
- `n::ClassicalNoise`: The classical noise.

# Returns
A function that takes a frequency `ω` and returns the corresponding spectrum value.
"""
spectrum(n::ClassicalNoise) = ω -> iszero(ω) ? zero(ω) : 2*n.T/ω

# Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
# statement makes sure the function does not divide by zero in case of `ω==0`.
"""
    spectrum(n::QuantumNoise)

Construct the spectrum of [`QuantumNoise`](@ref) at a given frequency.

# Arguments
- `n::QuantumNoise`: The quantum noise.

# Returns
A function that takes a frequency `ω` and returns the corresponding spectrum value.
"""
spectrum(n::QuantumNoise) = ω -> iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

# Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
# statement makes sure the function does not divide by zero in case of `ω==0`. It differs
# from `spectrum(n::QuantumNoise)` in the fact that the zero-point noise is removed.
"""
    spectrum(n::NoZeroQuantumNoise)

Construct the spectrum of [`NoZeroQuantumNoise`](@ref) at a given frequency.

# Arguments
- `n::NoZeroQuantumNoise`: The no-zero-point-fluctuations quantum noise.

# Returns
A function that takes a frequency `ω` and returns the corresponding spectrum value.
"""
spectrum(n::NoZeroQuantumNoise) = ω -> iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2) - sign(ω))
