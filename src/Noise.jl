"""
    abstract type Noise

An abstract type used to represent different kinds of environment noise.

Any user-defined subtype of `Noise` should implement a method `spectrum(::NoiseSubType, ω)`
which returns the spectrum of the noise at the given frequency ω.
"""
abstract type Noise end

(n::Noise)(ω) = spectrum(n, ω)

"""
    spectrum(::Noise)

Return the spectrum of the given noise as a function of frequency.

# Arguments
- `n::Noise`: The noise type.

# Returns
A function that takes a frequency `ω` and returns the corresponding spectrum value.
"""
spectrum(n::NT) where NT <: Noise = ω -> spectrum(n, ω)

"""
    struct ClassicalNoise{TT<:Real} <: Noise

A subtype of `Noise` used to represent classical noise, with spectrum
```math
\\mathcal{S}(\\omega) = \\frac{2k_\\mathrm{B}T}{\\hbar\\omega}
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

spectrum(n::ClassicalNoise, ω) = iszero(ω) ? zero(ω) : 2*n.T/ω

# Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
# statement makes sure the function does not divide by zero in case of `ω==0`.

spectrum(n::QuantumNoise, ω) = iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

# Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
# statement makes sure the function does not divide by zero in case of `ω==0`. It differs
# from `spectrum(n::QuantumNoise)` in the fact that the zero-point noise is removed.

spectrum(n::NoZeroQuantumNoise, ω) = iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2) - sign(ω))

"""
    function psd(J::AbstractSD, noise::Noise)

Calculate the Power Spectral Density (PSD) for a given environment spectral
density `J` and noise model `noise`.

# Arguments
- `J::AbstractSD`: The environment spectral density.
- `noise::Noise`: The noise model for the environment.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

# Returns
A function `psd(ω)` that calculates the PSD at a given frequency `ω`.
"""
function psd(J::AbstractSD, noise::Noise)
    psd(ω) = imag_memory_kernel_ft(J, ω)*noise(ω)
    return psd
end

psd(J::LorentzianSD, noise::ClassicalNoise) = ω -> 2*noise.T*π*sdoverω(J, ω)

psd(J::LorentzianSD, noise::QuantumNoise) = ω -> iszero(ω) ? 2*noise.T*π*sdoverω(J, 0) : π*sd(J, ω)*coth(ω/2/noise.T)

psd(J::LorentzianSD, noise::NoZeroQuantumNoise) = ω -> iszero(ω) ? 2*noise.T*π*sdoverω(J, 0) : π*sd(J, ω)*(coth(ω/2/noise.T) - sign(ω))