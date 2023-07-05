"""
```Julia
kernel(J::LorentzianSD)
```

Returns the specific damping kernel for a Lorentzian spectral density defined by the
parameters in `J`. The returned function depends on `ω`.
"""
kernel(J::LorentzianSD) = ω -> J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)

"""
```Julia
imagkernel(J::AbstractSD)
```

Returns the imaginary part of the damping kernel for a Generic spectral density real and anti-symmetric in ω.
Of this kind, we find Lorentzian spectral densities, and Polynomial spectral densities with ``s\\in\\mathbb{N}``.
The spectral density is defined by the parameters in `J`. The returned function depends on `ω`.
"""
imagkernel(J::AbstractSD) = ω -> π*s*J(ω)

"""
```Julia
psd(J::GenericSD, noise::Noise)
```

Returns the power spectral density depending on parameters `J` and `noise`. The returned
function depends on `ω`.
"""
function psd(J::AbstractSD, noise::Noise)
    imagK = imagkernel(J)
    n = spectrum(noise)
    psd(ω) = imagK(ω)*n(ω)
    return psd
end

"""
```Julia
psd(J::LorentzianSD, noise::ClassicalNoise)
```

Returns the analytical expression for power spectrum depending on Lorentzian spectral
density and Classical noise. The returned function depends on `ω`.
"""
psd(J::LorentzianSD, noise::ClassicalNoise) = ω -> 2*J.α*J.Γ*noise.T/((J.ω0^2 - ω^2)^2 + (J.Γ*ω)^2)