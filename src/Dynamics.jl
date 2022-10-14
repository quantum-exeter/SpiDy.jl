"""
```Julia
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[])
```

Returns `[sol.t, s, sinterp]`, that is,
- the vector `sol.t` of time steps at which the solutions are evaluated,
- the 3-vector solution `s[:, 1+(i-1)*3]`, `s[:, 2+(i-1)*3]`, `s[:, 3+(i-1)*3]` relative to the spin `i` is evaluated at times `sol.t`,
- the 3 functions `sinterp(t)[1+(i-1)*3]`, `sinterp(t)[2+(i-1)*3]`, `sinterp(t)[3+(i-1)*3]` are the interpolations of the relative solutions `s` found in the given time span.

The differential equation solver is built to account for Lorentzian spectral density.

Keyword arguments:
- `JH` is the Heisenberg coupling matrix. Note that this have to be a symmetric matrix with zero diagonal. The preset value is the additive identity of the UniformScaling type, `JH=zero(I)`.
- `S0` spin length set at default value 1/2, `S0=1/2`.
- `Bext` external magnetic field set as unit-vector along the z-axis as default, `Bext = [0, 0, 1]`.
- `saveat` is an option of the function `solve()` which allows to only save the solution at the points needed to evaluate the steady-state, i.e. at late times. Used to optimize memory management and speed of the solver. Default value is an empty list, `saveat=[]`, resulting in the solution being saved at optimal time steps within the time span.

# Examples
```julia-repl
julia> diffeqsolver(s0, tspan, J, bfields, matrix; saveat=((N*4÷5):1:N)*Δt)
```
"""
function diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[])

    N = div(length(s0), 3)
    u0 = vcat(s0, [0, 0, 0, 0, 0, 0])
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    Cω2v = zeros(3)
    
    function f(du, u, par, t)
        s = @view u[1:3*N] # @view does not allocate values. No hard copy, just reference.
        v = @view u[1+3*N:3+3*N]
        w = @view u[4+3*N:6+3*N]
        Beff = Bext + bn(t) + mul!(Cω2v, Cω2, v)
        for i in 1:N
            du[1+(i-1)*3:3+(i-1)*3] = -cross(s[1+(i-1)*3:3+(i-1)*3], Beff + sum([JH[i,j] * s[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N]))
        end
        du[1+3*N:3+3*N] = w
        du[4+3*N:6+3*N] = -(J.ω0^2)*v -J.Γ*w -J.α*sum([s[(1+(j-1)*3):(3+(j-1)*3)] for j in 1:N])
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, maxiters=Int(1e7), save_idxs=1:3*N, saveat=saveat)
    return sol
end

"""
```Julia
diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])
```

Returns `[sol.t, x, p, xinterp, pinterp]`, that is,
- the vector `sol.t` of time steps at which the solutions are evaluated,
- the vectors of the solutions positions `x` and momenta `p` evaluated at times `sol.t`,
- the functions xinterp(t), pinterp(t) interpolations of the solutions found in the given time span.

The differential equation solver is built to account for Lorentzian spectral density.

Keyword arguments:
- `Ω` harmonic oscillator bare frequency set as default to `Ω=1.0`.
- `saveat` is an option of the function `solve()` which allows to only save the solution at the points needed to evaluate the steady-state, i.e. at late times. Used to optimize memory management and speed of the solver. Default value is an empty list, `saveat=[]`, resulting in the solution being saved at optimal time steps within the time span.

# Examples
```julia-repl
julia> diffeqsolver(x0, p0, tspan, J, bfields, matrix; saveat=((N*4÷5):1:N)*Δt)
```
"""
function diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])

    u0 = [x0[1], x0[2], x0[3], p0[1], p0[2], p0[3], 0, 0, 0, 0, 0, 0]
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    Cω2v = zeros(3)
    
    function f(du, u, par, t)
        x = @view u[1:3] # @view does not allocate values. No hard copy, just reference.
        p = @view u[4:6]
        v = @view u[7:9]
        w = @view u[10:12]
        Beff = bn(t) + mul!(Cω2v, Cω2, v)
        du[1:3] = p
        du[4:6] = -(Ω^2)*x + Beff
        du[7:9] = w
        du[10:12] = -(J.ω0^2)*v -J.Γ*w + J.α*x
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, maxiters=Int(1e7), saveat=saveat)

    x = zeros(length(sol.t), 3)
    p = zeros(length(sol.t), 3)
    for n in 1:length(sol.t)
        x[n, :] = sol.u[n][1:3]
        p[n, :] = sol.u[n][4:6]
    end
    xinterp = t -> sol(t)[1:3]
    pinterp = t -> sol(t)[4:6]
    
    return sol.t, x, p, xinterp, pinterp
end
