"""
```Julia
GenericSD
```

Definition of the abstract type `GenericSD`.
"""
abstract type GenericSD end

"""
```Julia
LorentzianSD{T<:Real}
```

Returns a `LorentzianSD` structure of type `GenericSD` built by passing 3 `Real` values.
The values are ordered as `LorentzianSD(α, ω0, Γ)`.

# Examples
```julia-repl
julia> LorentzianSD(1., 3., 8.)
```
"""
struct LorentzianSD{T<:Real} <: GenericSD
    α::T
    ω0::T
    Γ::T
end

"""
```Julia
OhmicSD{T<:Real}
```

Returns a `OhmicSD` structure of type `GenericSD` built by passing 3 `Real` values.
The values are ordered as `OhmicSD(s, α, ωcut)`.

# Examples
```julia-repl
julia> OhmicSD(1., 3., 100.)
```
"""
struct PolySD{T<:Real} <: GenericSD
    s::T
    α::T
    ωcut::T
end

"""
```Julia
sd(J::GenericSD)
```

Defines the spectral density for generic shapes `GenericSD`.
"""
sd(J::GenericSD) = ω -> sdoverω(J)(ω)*ω

"""
```Julia
sdoverω(J::GenericSD)
```

Returns the spectral density divided by `ω` for generic shapes `GenericSD`.
"""
sdoverω(J::GenericSD) = ω -> sd(J)(ω)/ω

"""
```Julia
reorgenergy(J::GenericSD)
```

Returns the reorganization energy numerically integrated as ``\\int_0^\\infty \\text{sdoverω}(\\omega)d\\omega``.
"""
reorgenergy(J::GenericSD) = quadgk(sdoverω(J), 0.0, Inf)[1]

"""
```Julia
sdoverω(J::LorentzianSD)
```

Returns the spectral density divided by `ω` for `LorentzianSD` shapes which naturally defines `sd(J::LorentzianSD)`.
"""
sdoverω(J::LorentzianSD) = ω -> (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)

"""
```Julia
sdoverω(J::PolySD)
```

Returns the spectral density divided by `ω` for `PolySD` shapes which naturally defines `sdoverω(J::PolySD)`.
"""
sd(J::PolySD) = ω -> 2 * J.α * J.ωcut^(1-s) * ω^s * exp(-ω/J.ωcut)

"""
```Julia
reorgenergy(J::LorentzianSD)
```

Returns the analytical reorganization energy for `LorentzianSD` shapes.
"""
reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2


"""
```Julia
kernel(J::LorentzianSD)
```

Returns the specific damping kernel for a Lorentzian spectral density defined by the
parameters in `J`. It returns a function depending on `ω`.
"""
kernel(J::LorentzianSD) = ω -> J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)

"""
```Julia
psd(J::GenericSD, noise::Noise)
```

Returns the power spectral density depending on parameters `J` and `noise`. The returned
function depends on `ω`.
"""
function psd(J::GenericSD, noise::Noise)
    K = kernel(J)
    n = spectrum(noise)
    psd(ω) = imag(K(ω))*n(ω)
    return psd
end

"""
```Julia
psd(J::LorentzianSD, noise::ClassicalNoise)
```

Returns the analytical expression for power spectrum depending on Lorentzian spectral
density and Classical noise.
"""
psd(J::LorentzianSD, noise::ClassicalNoise) = ω -> 2*J.α*J.Γ*noise.T/((J.ω0^2 - ω^2)^2 + (J.Γ*ω)^2)