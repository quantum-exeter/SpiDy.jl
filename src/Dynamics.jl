"""
    function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], save_fields=false, projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves the dynamics of a system of ineracting spins under the influence of *both* local (unique to each spin)
and global (shared by all spins) stochastic noise from the environment.

# Arguments
- `s0`: Array of length `3N` specifying the initial conditions of the `N` spins. The order the initial consitions is first the `Sx,Sy,Sz` for the first spin, then for the second, and so on.
- `tspan`: The time span to solve the equations over, specified as a tuple `(tstart, tend)`.
- `Jlist::Vector{LorentzianSD}`: The list of spectral densities of the environment(s).
- `bfield`: An array of tuples of functions `Array{Tuple{Function, Function, Function}}` representing the time series of the noise for each environment.
- `bcoupling::Vector{Vector{T}} T <: Real`: A vector (of length the number of baths) of vectors (of length the number of spins), specifying if each bath couples to each spin.
- `matrix::Vector{TT} TT <: Coupling`: A vector of `Coupling` structs specifying the spin-environment coupling matrix for each environment.
- `JH=zero(I)`: (Optional) The spin-spin coupling matrix. Default is zero matrix (i.e. non-interacting spins).
- `S0=1/2`: (Optional) The spin quantum number. Default is 1/2.
- `Bext=[0, 0, 1]`: (Optional) The external magnetic field vector. Default is `[0, 0, 1]` (normalised length pointing in the `z` direction).
- `saveat=[]`: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.
- `save_fields=false`: (Optional) If true, also return the auxiliary fields encoding the environment memory.
- `projection=false`: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is `false`.
- `alg=Vern6()`: (Optional) The differential equation solver algorithm. Default is `Vern6()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-6`: (Optional) The absolute tolerance for the solver. Default is `1e-6`.
- `rtol=1e-6`: (Optional) The relative tolerance for the solver. Default is `1e-6`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

Note: Additional keyword arguments will be passed on to the ODE solver (see the `DifferentialEquations.jl` docs)

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(
    s0,
    tspan,
    Jlist::Vector{LorentzianSD},
    bfield,
    bcoupling::Vector{<:AbstractArray{T,1}} where {T<:Real},
    matrix::Vector{TT} where {TT<:Coupling};
    JH=zero(I),
    S0=1/2,
    Bext=[0, 0, 1],
    saveat=[],
    save_fields=false,
    projection=false,
    alg=Vern6(),
    atol=1e-6,
    rtol=1e-6,
    kwargs...
)
    N = div(length(s0), 3)
    if length(Jlist) != length(matrix) || length(Jlist) != length(bcoupling)
        throw(DimensionMismatch("The dimension of Jlist, bcoupling, and matrix must match."))
    end
    M = length(Jlist)
    u0 = [s0; zeros(6*M)]
    invsqrtS0 = 1/sqrt(S0)
    Cω = [matrix[i].C for i in 1:M]
    Cω2 = [matrix[i].C*transpose(matrix[i].C) for i in 1:M]
 
    b, Cb, Cω2v, Beff = dualcache(zeros(3)), dualcache(zeros(3)), dualcache(zeros(3)), dualcache(zeros(N, 3))
    params = (N, M, invsqrtS0, Bext, JH, Jlist, Cω, Cω2, bfield, bcoupling, b, Cb, Cω2v, Beff)
    prob = ODEProblem(_spin_time_step!, u0, tspan, params)
    condition(u, t, integrator) = true
    function affect!(integrator) # projection
        for n in 1:N
            integrator.u[1+(n-1)*3:3+(n-1)*3] ./= norm(integrator.u[1+(n-1)*3:3+(n-1)*3])
        end
        u_modified!(integrator, false)
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(false,false))
    skwargs = projection ? (callback=cb,) : NamedTuple()
    if save_fields
        save_idxs = 1:(3*N+6*M)
    else
        save_idxs = 1:3*N
    end
    sol = solve(prob, alg; abstol=atol, reltol=rtol, maxiters=Int(1e9), save_idxs=save_idxs, saveat=saveat, kwargs..., skwargs...)
    return sol
end

function _spin_time_step!(
    du,
    u,
    (N, M, invsqrtS0, Bext, JH, Jlist, Cω, Cω2, bfields, bcoupling, b, Cb, Cω2v, Beff),
    t
)
    b = get_tmp(b, u)
    Cb = get_tmp(Cb, u)
    Cω2v = get_tmp(Cω2v, u)
    Beff = get_tmp(Beff, u)

    s = @view u[1:3*N]
    v = @view u[1+3*N:3*N+3*M]
    w = @view u[1+3*N+3*M:3*N+6*M]
    ds = @view du[1:3*N]
    dv = @view du[1+3*N:3*N+3*M]
    dw = @view du[1+3*N+3*M:3*N+6*M]

    for i in 1:N
        Beff[i, :] .= Bext
    end

    for j in 1:M
        vj = @view v[1+(j-1)*3:3+(j-1)*3]
 
        for k in 1:3
            b[k] = bfields[j][k](t)
        end

        mul!(Cb, Cω[j], b)
        lmul!(invsqrtS0, Cb)
        mul!(Cω2v, Cω2[j], vj)
        Cb .+= Cω2v

        for i in 1:N
            for k in 1:3
                Beff[i, k] += bcoupling[j][i]*Cb[k]
            end
        end
    end

    for i in 1:N
        ds[1+(i-1)*3] = -(s[2+(i-1)*3]*Beff[i,3]-s[3+(i-1)*3]*Beff[i,2])
        ds[2+(i-1)*3] = -(s[3+(i-1)*3]*Beff[i,1]-s[1+(i-1)*3]*Beff[i,3])
        ds[3+(i-1)*3] = -(s[1+(i-1)*3]*Beff[i,2]-s[2+(i-1)*3]*Beff[i,1])

        for j in 1:N
            ds[1+(i-1)*3] += -(s[2+(i-1)*3]*JH[i,j]*s[3+(j-1)*3]-s[3+(i-1)*3]*JH[i,j]*s[2+(j-1)*3])
            ds[2+(i-1)*3] += -(s[3+(i-1)*3]*JH[i,j]*s[1+(j-1)*3]-s[1+(i-1)*3]*JH[i,j]*s[3+(j-1)*3])
            ds[3+(i-1)*3] += -(s[1+(i-1)*3]*JH[i,j]*s[2+(j-1)*3]-s[2+(i-1)*3]*JH[i,j]*s[1+(j-1)*3])
        end
    end

    dv .= w

    for i in 1:M
        vi = @view v[1+(i-1)*3:3+(i-1)*3]
        wi = @view w[1+(i-1)*3:3+(i-1)*3]

        for k in 1:3
            dw[k+3*(i-1)] = -(Jlist[i].ω0^2)*vi[k] -Jlist[i].Γ*wi[k]
        end

        for j in 1:N
            sj = @view s[1+3*(j-1):3+3*(j-1)]

            for k in 1:3
                dw[k+3*(i-1)] += -Jlist[i].α*bcoupling[i][j]*sj[k]
            end
        end
    end
end

"""
    function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], save_fields=false, projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves the dynamics of a system of ineracting spins under the influence of *either* local (unique to each spin)
or global (shared by all spins) stochastic noise from the environment.

# Arguments
- `s0`: Array of length `3N` specifying the initial conditions of the `N` spins. The order the initial consitions is first the `Sx,Sy,Sz` for the first spin, then for the second, and so on.
- `tspan`: The time span to solve the equations over, specified as a tuple `(tstart, tend)`.
- `J::LorentzianSD`: The spectral density of the noise acting on the spins (either local or shared depending on the value of `bfields`).
- `bfields`: For *local* baths, an array of tuples of functions `Array{Tuple{Function, Function, Function}}` representing the time series of the local stochastic field for each spin. For a *global* bath, a tuple of functions `Tuple{Function, Function, Function}` representing the time series of the global stochastic field shared by all the spins.
- `matrix::Coupling`: The spin-environment coupling matrix.
- `JH=zero(I)`: (Optional) The spin-spin coupling matrix. Default is zero matrix (i.e. non-interacting spins).
- `S0=1/2`: (Optional) The spin quantum number. Default is 1/2.
- `Bext=[0, 0, 1]`: (Optional) The external magnetic field vector. Default is `[0, 0, 1]` (normalised length pointing in the `z` direction).
- `saveat=[]`: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.
- `save_fields=false`: (Optional) If true, also return the auxiliary fields encoding the environment memory.
- `projection=false`: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is `false`.
- `alg=Vern6()`: (Optional) The differential equation solver algorithm. Default is `Vern6()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-6`: (Optional) The absolute tolerance for the solver. Default is `1e-6`.
- `rtol=1e-6`: (Optional) The relative tolerance for the solver. Default is `1e-6`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

Note: Additional keyword arguments will be passed on to the ODE solver (see the `DifferentialEquations.jl` docs)

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(
    s0,
    tspan,
    J::LorentzianSD,
    bfield,
    matrix::Coupling;
    JH=zero(I),
    S0=1/2,
    Bext=[0, 0, 1],
    saveat=[],
    save_fields=false,
    projection=false,
    alg=Vern6(),
    atol=1e-6,
    rtol=1e-6,
    kwargs...
)
    N = div(length(s0), 3)
    if length(bfield) == N && length(bfield[1]) == 3 # only local baths
        Jlist = repeat([J], N)
        bcoupling = [I(N)[i,:] for i in 1:N]
        matrix = repeat([matrix], N)
        return diffeqsolver(s0, tspan, Jlist, bfield, bcoupling, matrix; JH=JH, S0=S0, Bext=Bext, saveat=saveat, save_fields=save_fields, projection=projection, alg=alg, atol=atol, rtol=rtol, kwargs...)
    else # only shared bath
        return diffeqsolver(s0, tspan, [J], [bfield], [ones(N)], [matrix]; JH=JH, S0=S0, Bext=Bext, saveat=saveat, save_fields=save_fields, projection=projection, alg=alg, atol=atol, rtol=rtol, kwargs...)
    end
end

"""
    function diffeqsolver(x0, p0, tspan, Jlist::Vector{LorentzianSD}, bfield, bcoupling::Vector{<:AbstractArray{T,1}} where {T<:Real}, matrix::Vector{TT} where {TT<:Coupling}; JH=zero(I),  Ω=1.0, counter_term=true, saveat=[], save_fields=false, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves the dynamics of a system of ineracting oscillators under the influence of *both* local (unique to each oscillator)
and global (shared by all oscillators) stochastic noise from the environment.

# Arguments
- `x0`: Array of length `N` specifying the initial position of the `N` oscillators.
- `p0`: Array of length `N` specifying the initial momentum of the `N` oscillators.
- `tspan`: The time span to solve the equations over, specified as a tuple `(tstart, tend)`.
- `Jlist::Vector{LorentzianSD}`: The list of spectral densities of the environment(s).
- `bfield`: An array of tuples of functions `Array{Tuple{Function, Function, Function}}` representing the time series of the noise for each environment.
- `bcoupling::Vector{<:AbstractArray{T,1}} T <: Real`: A vector (of length the number of baths) of vectors (of length the number of oscillators), specifying if each bath couples to each oscillator.
- `matrix::Vector{TT} TT <: Coupling`: A vector of `Coupling` structs specifying the oscillator-environment coupling matrix for each environment.
- `JH=zero(I)`: (Optional) The oscillator-oscillator coupling matrix. Default is zero matrix (i.e. non-interacting oscillators).
- `Ω=1.0`: (Optional) The natural angular frequency of the harmonic oscillators (currently the same for all). Default is 1.
- `counter_term=true`: (Optional) Whether to include the counter-term or not. Default is true.
- `saveat=[]`: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.
- `save_fields=false`: (Optional) If true, also return the auxiliary fields encoding the environment memory.
- `alg=Vern6()`: (Optional) The differential equation solver algorithm. Default is `Vern6()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-6`: (Optional) The absolute tolerance for the solver. Default is `1e-6`.
- `rtol=1e-6`: (Optional) The relative tolerance for the solver. Default is `1e-6`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

Note: Additional keyword arguments will be passed on to the ODE solver (see the `DifferentialEquations.jl` docs)

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""

function diffeqsolver(
    x0,
    p0,
    tspan,
    Jlist::Vector{LorentzianSD},
    bfield,
    bcoupling::Vector{<:AbstractArray{T,1}} where {T<:Real},
    matrix::Vector{TT} where {TT<:Coupling};
    JH=zero(I),
    Ω=1.0,
    counter_term=true,
    saveat=[],
    save_fields=false,
    alg=Vern6(),
    atol=1e-6,
    rtol=1e-6,
    kwargs...
)
    N = div(length(x0), 3)
    if length(Jlist) != length(matrix) || length(Jlist) != length(bcoupling)
        throw(DimensionMismatch("The dimension of Jlist, bcoupling, and matrix must match."))
    end
    M = length(Jlist)
    u0 = [x0; p0; zeros(6*M)]
    Cω2 = []
    b = []
    Ωeff2 = repeat([Ω^2], N)
    for i in 1:M
        push!(Cω2, matrix[i].C*transpose(matrix[i].C))
        push!(b, t -> matrix[i].C*[bfield[i][1](t), bfield[i][2](t), bfield[i][3](t)]);
        if counter_term
            for j in 1:N
                Ωeff2[j] += 2*reorganisation_energy(Jlist[i])*bcoupling[i][j]
            end
        end
    end
    function f(du, u, (Cω2v, Beff), t)
        Cω2v = get_tmp(Cω2v, u)
        Beff = get_tmp(Beff, u)
        x = @view u[1:3*N]
        p = @view u[1+3*N:6*N]
        v = @view u[1+6*N:6*N+3*M]
        w = @view u[1+6*N+3*M:6*N+6*M]
        dx = @view du[1:3*N]
        dp = @view du[1+3*N:6*N]
        dv = @view du[1+6*N:6*N+3*M]
        dw = @view du[1+6*N+3*M:6*N+6*M]
        for i in 1:N
            Beff[i, :] .= 0
            for j in 1:M
                Beff[i, :] .+= bcoupling[j][i]*(b[j](t) + mul!(Cω2v, Cω2[j], v[1+(j-1)*3:3+(j-1)*3]))
            end
        end
        for i in 1:N
            dx[1+(i-1)*3:3+(i-1)*3] = p[1+(i-1)*3:3+(i-1)*3]
            dp[1+(i-1)*3:3+(i-1)*3] = -Ωeff2[i]*x[1+(i-1)*3:3+(i-1)*3] + Beff[i, :]
            for j in 1:N
                dp[1+(i-1)*3:3+(i-1)*3]  += JH[i,j] * x[(1+(j-1)*3):(3+(j-1)*3)]
            end
        end
        dv .= w
        for i in 1:M
            dw[1+3*(i-1):3+3*(i-1)] = -(Jlist[i].ω0^2)*v[1+3*(i-1):3+3*(i-1)] -Jlist[i].Γ*w[1+3*(i-1):3+3*(i-1)]
            for j in 1:N
                dw[1+3*(i-1):3+3*(i-1)] += Jlist[i].α*bcoupling[i][j]*x[1+3*(j-1):3+3*(j-1)]
            end
        end
    end
    if save_fields
        save_idxs = 1:(6*N+6*M)
    else
        save_idxs = 1:6*N
    end
    prob = ODEProblem(f, u0, tspan, (dualcache(zeros(3)), dualcache(zeros(N,3))))
    sol = solve(prob, alg; abstol=atol, reltol=rtol, maxiters=Int(1e9), save_idxs=save_idxs, saveat=saveat, kwargs...)
    return sol
end

"""
    function diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfield, matrix::Coupling; JH=zero(I), Ω=1.0, counter_term=true, saveat=[], save_fields=false, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves the dynamics of a system of ineracting oscillators under the influence of *either* local (unique to each oscillator)
or global (shared by all oscillators) stochastic noise from the environment.

# Arguments
- `x0`: Array of length `N` specifying the initial position of the `N` oscillators.
- `p0`: Array of length `N` specifying the initial momentum of the `N` oscillators.
- `tspan`: The time span to solve the equations over, specified as a tuple `(tstart, tend)`.
- `J::LorentzianSD`: The spectral density of the noise acting globally on the harmonic oscillators.
- `bfield`: A tuple of functions `Tuple{Function, Function, Function}` representing the time series of the global stochastic field shared by all the harmonic oscillators.
- `matrix::Coupling`: The harmonic oscillators-environment coupling matrix.
- `JH=zero(I)`: (Optional) The oscillator-oscillator coupling matrix. Default is zero matrix (i.e. non-interacting oscillators).
- `Ω=1.0`: (Optional) The natural angular frequency of the harmonic oscillators (currently the same for all). Default is 1.
- `counter_term=true`: (Optional) Whether to include the counter-term or not. Default is true.
- `saveat=[]`: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.
- `save_fields=false`: (Optional) If true, also return the auxiliary fields encoding the environment memory.
- `alg=Vern6()`: (Optional) The differential equation solver algorithm. Default is `Vern6()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-6`: (Optional) The absolute tolerance for the solver. Default is `1e-6`.
- `rtol=1e-6`: (Optional) The relative tolerance for the solver. Default is `1e-6`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

Note: Additional keyword arguments will be passed on to the ODE solver (see the `DifferentialEquations.jl` docs)

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(
    x0,
    p0,
    tspan,
    J::LorentzianSD,
    bfield,
    matrix::Coupling;
    JH=zero(I),
    Ω=1.0,
    counter_term=true,
    saveat=[],
    save_fields=false,
    alg=Vern6(),
    atol=1e-6,
    rtol=1e-6,
    kwargs...
)
    N = div(length(x0), 3)
    if length(bfield) == N && length(bfield[1]) == 3 # only local baths
        Jlist = repeat([J], N)
        bcoupling = [I(N)[i,:] for i in 1:N]
        matrix = repeat([matrix], N)
        return diffeqsolver(x0, p0, tspan, Jlist, bfield, bcoupling, matrix; JH=JH, Ω=Ω, counter_term=counter_term, saveat=saveat, save_fields=save_fields, alg=alg, atol=atol, rtol=rtol, kwargs...)
    else # only shared bath
        return diffeqsolver(x0, p0, tspan, [J], [bfield], [ones(N)], [matrix]; JH=JH, Ω=Ω, counter_term=counter_term, saveat=saveat, save_fields=save_fields, alg=alg, atol=atol, rtol=rtol, kwargs...)
    end
end