"""
diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves a system of coupled ordinary differential equations (ODEs) using numerical integration. This function supports **multiple baths in input**, local and shared.

Parameters:
- s0 : Array
    Initial conditions for the ODEs.
- tspan : Tuple or AbstractVector
    Time span for integration, given as a tuple (tstart, tend) or an array of time points.
- J : LorentzianSD
    Physical parameters for the J-coupling interaction.
- Jshared : LorentzianSD
    Physical parameters for the shared J-coupling interaction.
- bfields : Array{Tuple{Function, Function, Function}}
    Local stochastic field components as functions of time.
- bfieldshared : Tuple{Function, Function, Function}
    Shared stochastic field components as functions of time.
- matrix : Coupling
    Coupling matrix between spins.

Keyword Arguments:
- JH : Matrix, optional (default: zero(I))
    Matrix representing the J-coupling interaction strengths.
- S0 : Number, optional (default: 1/2)
    Spin quantum number.
- Bext : Array{Number}, optional (default: [0, 0, 1])
    External magnetic field components.
- saveat : AbstractVector, optional (default: [])
    Time points to save the solution.
- projection : Bool, optional (default: true)
    Flag indicating whether to perform projection on the solution.
- alg : Algorithm, optional (default: Tsit5())
    Integration algorithm to use.
- atol : Number, optional (default: 1e-3)
    Absolute tolerance for the integration.
- rtol : Number, optional (default: 1e-3)
    Relative tolerance for the integration.

Returns:
- sol : ODESolution
    Solution to the system of ODEs.
"""
function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(s0), 3)
    u0 = [s0; zeros(6*N+6)]
    Cω2 = matrix.C*transpose(matrix.C)
    b = []
    for i in 1:N
        push!(b, t -> matrix.C*[bfields[i][1](t), bfields[i][2](t), bfields[i][3](t)]);
    end
    bshared = t -> matrix.C*[bfieldshared[1](t), bfieldshared[2](t), bfieldshared[3](t)];
    function f(du, u, (Cω2v, Beff), t)
        Cω2v = get_tmp(Cω2v, u)
        Beff = get_tmp(Beff, u)
        s = @view u[1:3*N]
        v = @view u[1+3*N:6*N]
        w = @view u[1+6*N:9*N]
        vshared = @view u[1+9*N:3+9*N]
        wshared = @view u[4+9*N:6+9*N]
        for i in 1:N
            Beff[i, :] .= Bext + bshared(t) + b[i](t) + mul!(Cω2v, Cω2, vshared + v[1+(i-1)*3:3+(i-1)*3])
        end
        for i in 1:N
            du[1+(i-1)*3] = -(s[2+(i-1)*3]*Beff[i,3]-s[3+(i-1)*3]*Beff[i,2])
            du[2+(i-1)*3] = -(s[3+(i-1)*3]*Beff[i,1]-s[1+(i-1)*3]*Beff[i,3])
            du[3+(i-1)*3] = -(s[1+(i-1)*3]*Beff[i,2]-s[2+(i-1)*3]*Beff[i,1])
            for j in 1:N
                du[1+(i-1)*3] += -(s[2+(i-1)*3]*JH[i,j]*s[3+(j-1)*3]-s[3+(i-1)*3]*JH[i,j]*s[2+(j-1)*3])
                du[2+(i-1)*3] += -(s[3+(i-1)*3]*JH[i,j]*s[1+(j-1)*3]-s[1+(i-1)*3]*JH[i,j]*s[3+(j-1)*3])
                du[3+(i-1)*3] += -(s[1+(i-1)*3]*JH[i,j]*s[2+(j-1)*3]-s[2+(i-1)*3]*JH[i,j]*s[1+(j-1)*3])
            end
        end
        du[1+3*N:6*N] = w
        du[1+6*N:9*N] = -(J.ω0^2)*v -J.Γ*w -J.α*s
        du[1+9*N:3+9*N] = wshared
        du[4+9*N:6+9*N] = -(Jshared.ω0^2)*vshared -Jshared.Γ*wshared
        for i in 1:N
            du[4+9*N] += -Jshared.α*s[(1+(i-1)*3)]
            du[5+9*N] += -Jshared.α*s[(2+(i-1)*3)]
            du[6+9*N] += -Jshared.α*s[(3+(i-1)*3)]
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
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves a system of coupled ordinary differential equations (ODEs) using numerical integration. This function supports **a single bath input**, either local or shared.

For local baths:
- s0 : Array
    Initial conditions for the ODEs.
- tspan : Tuple or AbstractVector
    Time span for integration, given as a tuple (tstart, tend) or an array of time points.
- J : LorentzianSD
    Physical parameters for the J-coupling interaction.
- bfields : Array{Tuple{Function, Function, Function}}
    Local stochastic field components as functions of time.
- matrix : Coupling
    Coupling matrix between spins.

For shared bath:
- s0 : Array
    Initial conditions for the ODEs.
- tspan : Tuple or AbstractVector
    Time span for integration, given as a tuple (tstart, tend) or an array of time points.
- J : LorentzianSD
    Physical parameters for the J-coupling interaction.
- bfields : Tuple{Function, Function, Function}
    Shared stochastic field components as functions of time.
- matrix : Coupling
    Coupling matrix between spins.

Keyword Arguments:
- JH : Matrix, optional (default: zero(I))
    Matrix representing the J-coupling interaction strengths.
- S0 : Number, optional (default: 1/2)
    Spin quantum number.
- Bext : Array{Number}, optional (default: [0, 0, 1])
    External magnetic field components.
- saveat : AbstractVector, optional (default: [])
    Time points to save the solution.
- projection : Bool, optional (default: true)
    Flag indicating whether to perform projection on the solution.
- alg : Algorithm, optional (default: Tsit5())
    Integration algorithm to use.
- atol : Number, optional (default: 1e-3)
    Absolute tolerance for the integration.
- rtol : Number, optional (default: 1e-3)
    Relative tolerance for the integration.

Returns:
- sol : ODESolution
    Solution to the system of ODEs.
"""
function diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(s0), 3)
    if length(bfields) == N && length(bfields[1]) == 3 # only local baths
        Jshared = LorentzianSD(0, 0, 0)
        bfieldshared = [t -> 0, t -> 0, t -> 0]
        return diffeqsolver(s0, tspan, J, Jshared, bfields, bfieldshared, matrix; JH=JH, S0=S0, Bext=Bext, saveat=saveat, projection=projection, alg=alg, atol=atol, rtol=rtol)
    else # only shared bath
        Jlocal = LorentzianSD(0, 0, 0)
        bfieldslocal = []
        for _ in 1:N
            push!(bfieldslocal, [t -> 0, t -> 0, t -> 0])
        end
        return diffeqsolver(s0, tspan, Jlocal, J, bfieldslocal, bfields, matrix; JH=JH, S0=S0, Bext=Bext, saveat=saveat, projection=projection, alg=alg, atol=atol, rtol=rtol)
    end
end

"""
diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), Ω=1.0, saveat=[], alg=Tsit5(), atol=1e-3, rtol=1e-3)

Solves a system of coupled ordinary differential equations (ODEs) using numerical integration.

Parameters:
- x0 : Array
    Initial positions for the ODEs.
- p0 : Array
    Initial momenta for the ODEs.
- tspan : Tuple or AbstractVector
    Time span for integration, given as a tuple (tstart, tend) or an array of time points.
- J : LorentzianSD
    Physical parameters for the J-coupling interaction.
- bfields : Tuple{Function, Function, Function}
    Magnetic field components as functions of time.
- matrix : Coupling
    Coupling matrix.

Keyword Arguments:
- JH : Matrix, optional (default: zero(I))
    Matrix representing the J-coupling interaction strengths.
- Ω : Number, optional (default: 1.0)
    Angular frequency.
- saveat : AbstractVector, optional (default: [])
    Time points to save the solution.
- alg : Algorithm, optional (default: Tsit5())
    Integration algorithm to use.
- atol : Number, optional (default: 1e-3)
    Absolute tolerance for the integration.
- rtol : Number, optional (default: 1e-3)
    Relative tolerance for the integration.

Returns:
- sol : ODESolution
    Solution to the system of ODEs.
"""
function diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), Ω=1.0, saveat=[], alg=Tsit5(), atol=1e-3, rtol=1e-3)
    N = div(length(x0), 3)
    u0 = [x0; p0; [0, 0, 0, 0, 0, 0]]
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    Cω2v = zeros(3)
    Beff = zeros(3)
    function f(du, u, par, t)
        x = @view u[1:3*N]
        p = @view u[1+3*N:6*N]
        v = @view u[1+6*N:3+6*N]
        w = @view u[4+6*N:6+6*N]
        Beff .= bn(t) + mul!(Cω2v, Cω2, v)
        for i in 1:N
            du[1+(i-1)*3:3+(i-1)*3] = p[1+(i-1)*3:3+(i-1)*3]
            du[1+3*N+(i-1)*3:3+3*N+(i-1)*3] = -(Ω^2)*x[1+(i-1)*3:3+(i-1)*3] + Beff + sum([JH[i,j] * x[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N])
        end
        du[1+6*N:3+6*N] = w
        du[4+6*N:6+6*N] = -(J.ω0^2)*v -J.Γ*w +J.α*sum([x[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N])
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, alg, abstol=atol, reltol=rtol, maxiters=Int(1e7), save_idxs=1:6*N, saveat=saveat)
    return sol
end
