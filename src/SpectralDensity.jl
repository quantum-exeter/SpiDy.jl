## Spectral density type ##
abstract type SpectralDensity end

## Spectral density and spectral-density-divided-by-ω for generic shapes ##
sd(J::SpectralDensity) = ω -> sdoverω(J)(ω)*ω
sdoverω(J::SpectralDensity) = ω -> sd(J)(ω)/ω

## Reorganization energy numerically integrated ∫_0^inf{spectral density}dω ##
reorgenergy(J::SpectralDensity) = quadgk(sdoverω(J), 0.0, Inf)[1]

## Lorentzian Spectral Density ##
struct LorentzianSD{T<:Real} <: SpectralDensity
  α::T
  ω0::T
  Γ::T
end

## Spectral density divided by ω which naturally defines sd(J::SpectralDensity) ##
sdoverω(J::LorentzianSD) = ω -> (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)

## Analytical reorganization energy for Lorentzian spectral density ##
reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2

## Specific damping kernel for a Lorentzian spectral density ##
## Depends on ω and J ##
function damping_kernel_frequency(J::LorentzianSD)
  function K(ω)
    return J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)
  end
  return K
end