## Power spectrum depending on spectral density and noise ##
function PSD(J::SpectralDensity, noise::Noise)
  K = damping_kernel_frequency(J)
  n = spectrum(noise)
  psd(ω) = imag(K(ω))*n(ω)
  return psd
end