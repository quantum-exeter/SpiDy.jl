#### Spectral densities types ####
abstract type SpectralDensity end

sd(J::SpectralDensity) = ω -> sdinvω(J,ω)*ω
sd(J::SpectralDensity, ω) = sdinvω(J)(ω)

sdinvω(J::SpectralDensity) = ω -> sd(J,ω)/ω
sdinvω(J::SpectralDensity, ω) = sd(J)(ω)

reorgenergy(J::SpectralDensity) = quadgk(ω -> sdinvω(J,ω), 0.0, Inf)[1]

# function damping_kernel_time(J::SpectralDensity)
#   function K(τ)
#     if τ <= zero(τ)
#       return zero(τ)
#     else
#       return 2*quadgk(ω -> sd(J,ω)*sin(ω*τ), 0.0, Inf)[1]
#     end
#   end
#   return Κ
# end

# function damping_kernel_frequency()
# end

## Lorentzian Spectral Density ##
struct LorentzianSD{T<:Real} <: SpectralDensity
  α::T
  ω0::T
  Γ::T
end

sdinvω(J::LorentzianSD) = ω -> (J.α*J.Γ/π)/((ω^2 - J.ω0^2)^2 + (J.Γ*ω)^2)

reorgenergy(J::LorentzianSD) = (J.α/J.ω0^2)/2

function damping_kernel_time(J::LorentzianSD)
  function K(τ)
    if τ <= zero(τ)
      return zero(τ)
    else
      ω1 = sqrt(J.ω0^2 - (J.Γ/2)^2)
      return J.α*exp(-J.Γ*τ/2)*sin(ω1*τ)/ω1
    end
  end
  return K
end

function damping_kernel_frequency(J::LorentzianSD)
  function K(ω)
    return J.α/(J.ω0^2 - ω^2 - 1im*ω*J.Γ)
  end
  return K
end