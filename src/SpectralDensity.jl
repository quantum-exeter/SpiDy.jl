## Spectral density type ##
abstract type GenericSD end

## Lorentzian Spectral Density structure##
struct LorentzianSD{T<:Real} <: GenericSD
  α::T
  ω0::T
  Γ::T
end

"""
sd(J::GenericSD)

Spectral density for generic shapes.
"""
sd(J::GenericSD) = ω -> sdoverω(J)(ω)*ω

"""
sdoverω(J::GenericSD)

Spectral-density-divided-by-ω for generic shapes.
"""
sdoverω(J::GenericSD) = ω -> sd(J)(ω)/ω

"""
reorgenergy(J::GenericSD)

Reorganization energy numerically integrated as ``\\int_0^\\infty \\text{sdoverω}(\\omega)d\\omega``.
"""
reorgenergy(J::GenericSD) = quadgk(sdoverω(J), 0.0, Inf)[1]

"""
sdoverω(J::LorentzianSD)

Spectral density divided by ω which naturally defines sd(J::LorentzianSD).
"""
sdoverω(J::LorentzianSD) = ω -> (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)

"""
reorgenergy(J::LorentzianSD)

Analytical reorganization energy for Lorentzian spectral density.
"""
reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2


"""
kernel(J::LorentzianSD)

Specific damping kernel for a Lorentzian spectral density defined by the
parameters in J. It returns a function depending on ω.
"""
kernel(J::LorentzianSD) = ω -> J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)

"""
psd(J::GenericSD, noise::Noise)

Power spectral density depending on parameters J and noise. It
returns a function of ω.
"""
function psd(J::GenericSD, noise::Noise)
  K = kernel(J)
  n = spectrum(noise)
  psd(ω) = imag(K(ω))*n(ω)
  return psd
end

"""
psd(J::LorentzianSD, noise::ClassicalNoise)

Returns the analytical expression for power spectrum depending on Lorentzian spectral
density and Classical noise.
"""
psd(J::LorentzianSD, noise::ClassicalNoise) = ω -> 2*J.α*J.Γ*noise.T/((J.ω0^2 - ω^2)^2 + (J.Γ*ω)^2)