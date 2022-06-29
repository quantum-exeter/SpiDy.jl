"""
psd(J::SpectralDensity, noise::Noise)

Returns the power spectrum depending on spectral density and noise.
"""
function psd(J::SpectralDensity, noise::Noise)
  K = damping_kernel_frequency(J)
  n = spectrum(noise)
  psd(ω) = imag(K(ω))*n(ω)
  return psd
end

psd(J::LorentzianSD, noise::ClassicalNoise) = ω -> 2*J.α*J.Γ*noise.T/((J.ω0^2 - ω^2)^2 + (J.Γ*ω)^2)

"""
b_field(N, Δt, J::SpectralDensity, noise::Noise, distro=nothing)

Returns the stochastic field b(t). It is evaluated using the Lorentzian spectral
density defined by the parameters J, the classical/quantum/Barker noise. The sampling
of the stochastic noise is done in frequency space. The default stochastic noise
is white noise having Gaussian distribution but different distributions can be specified.
N defines the number of steps and Δt defines the time step.
"""
function b_field(N, Δt, J::SpectralDensity, noise::Noise, distro=nothing)
  if isnothing(distro)
    distro = Normal(0., 1/sqrt(Δt))
  end
  # distrosample = rand(distro, N)
  # distrosamplefft = rfft(distrosample)
  psdfft = PSD(J, noise).(2π*rfftfreq(N, 1/Δt)) # NB: rfftfreq(N, 1/Δt) takes the frequency step in input!
  distrosamplefft = (randn(length(psdfft)) + 1im*randn(length(psdfft)))/sqrt(2)
  bfft = psdfft.*distrosamplefft
  return irfft(bfft, N)
end