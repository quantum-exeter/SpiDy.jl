"""
bfield(N, Δt, J::GenericSD, noise::Noise, distro=nothing)

Returns the stochastic field b(t). It is evaluated using the Lorentzian spectral
density defined by the parameters J, the classical/quantum/quantum-no-zero-energy noise.
The sampling of the stochastic noise is done in frequency space. The default stochastic noise
is white noise having Gaussian distribution but different distributions can be specified.
N defines the number of steps and Δt defines the time step.
"""
function bfield(N, Δt, J::GenericSD, noise::Noise, distro=Normal(0., 1/sqrt(Δt)))
    distrosample = rand(distro, N)
    distrosamplefft = rfft(distrosample)
    psdfft = 2π*psd(J, noise).(2π*rfftfreq(N, 1/Δt)) # NB: rfftfreq(N, 1/Δt) takes the frequency step in input!
    bfft = sqrt.(psdfft).*distrosamplefft
    return irfft(bfft, N)
end