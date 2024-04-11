# Returns the stochastic field `b(t)`. It is evaluated using the Lorentzian spectral
# density defined by the parameters `J`, the classical/quantum/quantum-no-zero-energy noise.
# The sampling of the stochastic noise is done in frequency space. The default stochastic noise
# is white noise having Gaussian distribution but different distributions can be specified.
# `N` defines the number of steps and `Δt` defines the time step.
"""
    bfield(N, Δt, J::AbstractSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true)

Generate a stochastic field realisation time series (of length `N` and spacing `Δt`)
based on the given noise model `noise` and spectral density `J`.

# Arguments
- `N`: The number of time steps.
- `Δt`: The time step size.
- `J::AbstractSD`: The environment spectral density.
- `noise::Noise`: The noise model for the environment.
- `distro=Normal(0., 1/sqrt(Δt))`: (Optional) The distribution of noise samples. Default is a normal distribution with mean 0 and standard deviation `1/sqrt(Δt)`.
- `interpolation=true`: (Optional) Specifies whether to use linear interpolation for the stochastic field time series. Default is `true`.

Note: The [`AbstractSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.AbstractSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

# Returns
A time series of the stochastic field values.
"""
function bfield(N, Δt, J::AbstractSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true)
    distrosample = rand(distro, N)
    distrosamplefft = rfft(distrosample)
    psdfft = psd(J, noise).(2π*rfftfreq(N, 1/Δt)) # NB: rfftfreq(N, 1/Δt) takes the frequency step in input!
    bfft = sqrt.(psdfft).*distrosamplefft
    b = irfft(bfft, N)
    return interpolation ? linear_interpolation(LinRange(0, N*Δt, N), b) : b
end