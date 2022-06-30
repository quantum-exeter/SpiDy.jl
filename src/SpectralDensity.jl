## Spectral density type ##
abstract type GenericSD end

## Spectral density and spectral-density-divided-by-ω for generic shapes ##
sd(J::GenericSD) = ω -> sdoverω(J)(ω)*ω
sdoverω(J::GenericSD) = ω -> sd(J)(ω)/ω

## Reorganization energy numerically integrated ∫_0^inf{spectral density}dω ##
reorgenergy(J::GenericSD) = quadgk(sdoverω(J), 0.0, Inf)[1]

## Lorentzian Spectral Density ##
struct LorentzianSD{T<:Real} <: GenericSD
  α::T
  ω0::T
  Γ::T
end

## Spectral density divided by ω which naturally defines sd(J::GenericSD) ##
sdoverω(J::LorentzianSD) = ω -> (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)

## Analytical reorganization energy for Lorentzian spectral density ##
reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2

## Specific damping kernel for a Lorentzian spectral density ##
## Depends on ω and J ##
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