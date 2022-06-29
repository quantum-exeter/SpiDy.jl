"""
PSD(J::SpectralDensity, noise::Noise)

Returns the power spectrum depending on spectral density and noise.
"""
function PSD(J::SpectralDensity, noise::Noise)
  K = damping_kernel_frequency(J)
  n = spectrum(noise)
  psd(ω) = imag(K(ω))*n(ω)
  return psd
end

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
  distrosample = rand(distro, N)
  distrosamplefft = rfft(distrosample)
  freqsample = PSD(J, noise).(2π*rfftfreq(N, 1/Δt)) # the dot evaluates the function at all the array elements
  bfft = freqsample.*distrosamplefft
  return irfft(bfft, N)
end