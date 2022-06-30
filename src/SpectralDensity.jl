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

Reorganization energy numerically integrated as ``\\int_0^\\infty \\text{psd}(\\omega)d\\omega``.
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
damping_kernel_frequency(J::LorentzianSD)

Specific damping kernel for a Lorentzian spectral density. It returns
a function depending on ω and J.
"""
damping_kernel_frequency(J::LorentzianSD) = ω -> J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)

"""
psd(J::GenericSD, noise::Noise)

Returns the power spectrum depending on spectral density and noise.
"""
function psd(J::GenericSD, noise::Noise)
  K = damping_kernel_frequency(J)
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