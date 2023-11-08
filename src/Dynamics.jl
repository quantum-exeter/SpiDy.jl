"""
    function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

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
- `projection=true`: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is `true`.
- `alg=Tsit5()`: (Optional) The differential equation solver algorithm. Default is `Tsit5()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-3`: (Optional) The absolute tolerance for the solver. Default is `1e-3`.
- `rtol=1e-3`: (Optional) The relative tolerance for the solver. Default is `1e-3`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(s0, tspan, Jlist::Vector{LorentzianSD}, bfield, bcoupling::Vector{<:AbstractArray{T,1}} where {T<:Real}, matrix::Vector{TT} where {TT<:Coupling}; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(s0), 3)
    if length(Jlist) != length(matrix) || length(Jlist) != length(bcoupling)
        throw(DimensionMismatch("The dimension of Jlist, bcoupling, and matrix must match."))
    end
    M = length(Jlist)
    u0 = [s0; zeros(6*M)]
    Cω2 = []
    b = []
    for i in 1:M
        push!(Cω2, matrix[i].C*transpose(matrix[i].C))
        push!(b, t -> matrix[i].C*[bfield[i][1](t), bfield[i][2](t), bfield[i][3](t)]);
    end
    function f(du, u, (Cω2v, Beff), t)
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
            for j in 1:M
                Beff[i, :] .+= bcoupling[j][i]*(b[j](t) + mul!(Cω2v, Cω2[j], v[1+(j-1)*3:3+(j-1)*3]))
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
            dw[1+3*(i-1):3+3*(i-1)] = -(Jlist[i].ω0^2)*v[1+3*(i-1):3+3*(i-1)] -Jlist[i].Γ*w[1+3*(i-1):3+3*(i-1)]
            for j in 1:N
                dw[1+3*(i-1):3+3*(i-1)] += -Jlist[i].α*bcoupling[i][j]*s[1+3*(j-1):3+3*(j-1)]
            end
        end
    end
    prob = ODEProblem(f, u0, tspan, (dualcache(zeros(3)), dualcache(zeros(N,3))))
    condition(u, t, integrator) = true
    function affect!(integrator) # projection
        for n in 1:N
            integrator.u[1+(n-1)*3:3+(n-1)*3] ./= norm(integrator.u[1+(n-1)*3:3+(n-1)*3])
        end
        u_modified!(integrator, false)
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(false,false))
    skwargs = projection ? (callback=cb,) : NamedTuple()
    sol = solve(prob, alg, abstol=atol, reltol=rtol, maxiters=Int(1e7), save_idxs=1:3*N, saveat=saveat; skwargs...)
    return sol
end

"""
    function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

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
- `projection=true`: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is `true`.
- `alg=Tsit5()`: (Optional) The differential equation solver algorithm. Default is `Tsit5()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-3`: (Optional) The absolute tolerance for the solver. Default is `1e-3`.
- `rtol=1e-3`: (Optional) The relative tolerance for the solver. Default is `1e-3`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(s0, tspan, J::LorentzianSD, bfield, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(s0), 3)
    if length(bfield) == N && length(bfield[1]) == 3 # only local baths
        Jlist = repeat([J], N)
        bcoupling = [I(N)[i,:] for i in 1:N]
        matrix = repeat([matrix], N)
        return diffeqsolver(s0, tspan, Jlist, bfield, bcoupling, matrix; JH=JH, S0=S0, Bext=Bext, saveat=saveat, projection=projection, alg=alg, atol=atol, rtol=rtol)
    else # only shared bath
        return diffeqsolver(s0, tspan, [J], [bfield], [ones(N)], [matrix]; JH=JH, S0=S0, Bext=Bext, saveat=saveat, projection=projection, alg=alg, atol=atol, rtol=rtol)
    end
end

"""
    function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves the dynamics of a system of ineracting spins under the influence of *either* local (unique to each spin)
or global (shared by all spins) stochastic noise from the environment.

# Arguments
- `x0`: Array of length `N` specifying the initial position of the `N` oscillators.
- `p0`: Array of length `N` specifying the initial momentum of the `N` oscillators.
- `tspan`: The time span to solve the equations over, specified as a tuple `(tstart, tend)`.
- `J::LorentzianSD`: The spectral density of the noise acting globally on the harmonic oscillators.
- `bfields`: A tuple of functions `Tuple{Function, Function, Function}` representing the time series of the global stochastic field shared by all the harmonic oscillators.
- `matrix::Coupling`: The harmonic oscillators-environment coupling matrix.
- `JH=zero(I)`: (Optional) The oscillator-oscillator coupling matrix. Default is zero matrix (i.e. non-interacting oscillators).
- `Ω=1`: (Optional) The natural angular frequency of the harmonic oscillators (currently the same for all). Default is 1.
- `counter_term=true`: (Optional) Whether to include the counter-term or not. Default is true.
- `saveat=[]`: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.
- `projection=true`: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is `true`.
- `alg=Tsit5()`: (Optional) The differential equation solver algorithm. Default is `Tsit5()`. See the `DifferentialEquations.jl` docs for [choices](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
- `atol=1e-3`: (Optional) The absolute tolerance for the solver. Default is `1e-3`.
- `rtol=1e-3`: (Optional) The relative tolerance for the solver. Default is `1e-3`.

Note: The [`LorentzianSD`](https://quantum-exeter.github.io/SpectralDensities.jl/stable/reference/#SpectralDensities.LorentzianSD) type is provided by the [SpectralDensities.jl](https://github.com/quantum-exeter/SpectralDensities.jl) package.

# Returns
An [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/) struct from `DifferentialEquations.jl` containing the solution of the equations of motion.
"""
function diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), Ω=1.0, counter_term=true, saveat=[], alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(x0), 3)
    u0 = [x0; p0; [0, 0, 0, 0, 0, 0]]
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    Cω2v = zeros(3)
    Beff = zeros(3)
    Ωeff2 = counter_term ? Ω^2 + 2*reorganisation_energy(J) : Ω^2
    function f(du, u, par, t)
        x = @view u[1:3*N]
        p = @view u[1+3*N:6*N]
        v = @view u[1+6*N:3+6*N]
        w = @view u[4+6*N:6+6*N]
        Beff .= bn(t) + mul!(Cω2v, Cω2, v)
        for i in 1:N
            du[1+(i-1)*3:3+(i-1)*3] = p[1+(i-1)*3:3+(i-1)*3]
            du[1+3*N+(i-1)*3:3+3*N+(i-1)*3] = -Ωeff2*x[1+(i-1)*3:3+(i-1)*3] + Beff + sum([JH[i,j] * x[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N])
        end
        du[1+6*N:3+6*N] = w
        du[4+6*N:6+6*N] = -(J.ω0^2)*v -J.Γ*w +J.α*sum([x[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N])
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, alg, abstol=atol, reltol=rtol, maxiters=Int(1e7), save_idxs=1:6*N, saveat=saveat)
    return sol
end
