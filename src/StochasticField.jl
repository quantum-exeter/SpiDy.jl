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
function bfield(
    N,
    Δt,
    J::AbstractSD,
    noise::Noise;
    distro=Normal(0., 1/sqrt(Δt)),
    interpolation=true
)
    distrosample = rand(distro, N)
    distrosamplefft = rfft(distrosample)
    psdfft = psd(J, noise).(2π*rfftfreq(N, 1/Δt)) # NB: rfftfreq(N, 1/Δt) takes the frequency step in input!
    bfft = sqrt.(psdfft).*distrosamplefft
    b = irfft(bfft, N)
    return interpolation ? cubic_spline_interpolation(LinRange(0, N*Δt, N), b) : b
end

"""
    bfield(tend, J::AbstractSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true, atol=0, rtol=sqrt(eps()), maxiters=Int(1e6))

Generate a stochastic field realisation time series of length `tend` given the noise model
`noise` and spectral density `J`. The optional parameters `atol`, `rtol` control the
accuracy of the power spectral density of the generated noise.

# Arguments
- `tend`: The total evolution time.
- `J::AbstractSD`: The environment spectral density.
- `noise::Noise`: The noise model for the environment.
- `atol=0`: (Optional) The absolute tolerance for the spectral density estimation. Default is `0`.
- `rtol=sqrt(eps())`: (Optional) The relative tolerance for the spectral density estimation. Default is `sqrt(eps())`.
- `maxiters=Int(1e6)`: (Optional) The maximum number of iterations for the spectral density estimation. Default is `1e6`.
- `interpolation=true`: (Optional) Specifies whether to use linear interpolation for the stochastic field time series. Default is `true`.

Note: The parameters used for generating the stochastic field are estimated using the `estimate_bfield_parameters` function. See the documentation for more details.

# Returns
A time series of the stochastic field values.
"""
function bfield(
    tend,
    J::AbstractSD,
    noise::Noise;
    atol=0,
    rtol=1e-6,
    maxiters=Int(1e6),
    interpolation=true,
)
    dt, tmin = estimate_bfield_parameters(J, noise; atol=atol, rtol=rtol, maxiters=maxiters)
    tmax = max(tmin, tend)
    dt = iszero(dt) ? tmax/2 : dt
    N = 2*ceil(Int, tmax/dt)
    return bfield(N, dt, J, noise; interpolation=interpolation)
end

"""
    estimate_bfield_parameters(J::AbstractSD, noise::Noise; atol=0, rtol=sqrt(eps()), maxiters=Int(1e6))

Estimate the minimum time step size `Δt` and minimum total evolution time `tmin` required to
generate a stochastic field realisation time series based on the given noise model `noise`
and spectral density `J`, such that the power spectral density of the generated noise matches
the target power spectral density within absolute `atol` and relative tolerances `rtol`.

# Arguments
- `J::AbstractSD`: The environment spectral density.
- `noise::Noise`: The noise model for the environment.
- `atol=0`: (Optional) The absolute tolerance for the spectral density estimation. Default is `0`.
- `rtol=sqrt(eps())`: (Optional) The relative tolerance for the spectral density estimation. Default is `sqrt(eps())`.
- `maxiters=Int(1e6)`: (Optional) The maximum number of iterations for the spectral density estimation. Default is `1e6`.

# Returns
The minimum time step size `Δt` and minimum total evolution time `tmin`.
"""
function estimate_bfield_parameters(
    J::AbstractSD,
    noise::Noise;
    atol=0,
    rtol=1e-6,
    maxiters=Int(1e6)
)
    f = psd(J, noise)
    I, _, segs = quadgk_segbuf(f, 0.0, Inf; order=2)

    dω = _min_diff_seg(segs)
    ωcutoff = _min_seg_end(segs)

    # handle the case of the psd being exactly zero (e.g. classical T = 0 noise)
    if isinf(dω) || isinf(ωcutoff)
        return 0.0, 0.0
    end

    counter = 0
    while true
        counter > maxiters && throw(error("Failed to converge in $(maxiters) iterations"))
        It, _ = quadgk(f, 0.0, ωcutoff; order=2)
        Idiff = abs(It - I)
        if (Idiff < atol) || (Idiff < rtol*I)
            break
        end
        ωcutoff += 10*dω
        counter += 1
    end
    ωcutoff = 2*ωcutoff # ensure that the cutoff frequency is below the Nyquist frequency

    dt = 2π/ωcutoff
    tmin = 2π/dω

    return dt, tmin
end

_inv_coord(t) = t*inv(1 - t)

_min_seg_end(segs) = minimum([_inv_coord(s.b) for s in segs])

_min_diff_seg(segs) = minimum([abs(_inv_coord(s.b) - _inv_coord(s.a)) for s in segs])