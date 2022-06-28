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

function b_field(N, Δt, J::SpectralDensity, noise::Noise, distro=nothing)
  if isnothing(distro)
    distro = Normal(0., 1/sqrt(Δt))
  end
  distrosample = rand(distro, N)
  distrosamplefft = rfft(distrosample)
  freqsample = PSD(J, noise).(2π*rfftfreq(N, Δt)) # the dot evaluates the function at all the array elements
  bfft = freqsample.*distrosamplefft
  return irfft(bfft, N)
end