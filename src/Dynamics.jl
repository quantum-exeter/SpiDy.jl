"""
```Julia
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1], saveat=[])
```

Returns `[sol.t, s, sinterp]`, that is, the vector `sol.t` of time steps at which the solutions are evaluated,
the 3-vector of the solutions `s[1]`, `s[2]`, `s[3]` evaluated at times `sol.t`,
the 3 functions sinterp(t)[i] interpolations of the solutions `s[i]` found in the given time span.

The differential equation solver is built to account for Lorentzian spectral density.

Keyword arguments:
- `S0` spin length set at default value 1/2, `S0=1/2`.
- `Bext` external magnetic field set as unit-vector along the z-axis as default, `Bext = [0, 0, 1]`
- `saveat` is an option of the function `solve()` which allows to only save the solution at the points needed to evaluate the steady-state,
i.e. at late times. Used to optimize memory management and speed of the solver. Default value is an empty list, `saveat=[]`, resulting
in the solution saved at optimal time steps withing the entire time span.

# Examples
```julia-repl
julia> diffeqsolver(s0, tspan, J, bfields, matrix; saveat=((N*4÷5):1:N)*Δt)
```
"""
function diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1], saveat=[])

    u0 = [s0[1], s0[2], s0[3], 0, 0, 0, 0, 0, 0]
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    
    function f(du, u, par, t)
        s = @view u[1:3] # @view does not allocate values. No hard copy, just reference.
        v = @view u[4:6]
        w = @view u[7:9]
        du[1:3] = -cross(s, Bext + bn(t) + Cω2*v)
        du[4:6] = w
        du[7:9] = -(J.ω0^2)*v -J.Γ*w -J.α*s
    end
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Vern7(), abstol=1e-8, reltol=1e-8, maxiters=Int(1e7), saveat=saveat)

    s = zeros(length(sol.t), 3)
    for n in 1:length(sol.t)
        s[n, :] = sol.u[n][1:3]
    end
    sinterp = t -> sol(t)[1:3]
    
    return sol.t, s, sinterp
end

"""
```Julia
diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])
```

Returns `[sol.t, x, p, xinterp, pinterp]`, that is, the vector `sol.t` of time steps at which the solutions are evaluated,
the vectors of the solutions positions `x` and momenta `p` evaluated at times `sol.t`,
the functions xinterp(t), pinterp(t) interpolations of the solutions found in the given time span.

The differential equation solver is built to account for Lorentzian spectral density.

Keyword arguments:
- `Ω` harmonic oscillator bare frequency set as default to `Ω=1.0`
- `saveat` is an option of the function `solve()` which allows to only save the solution at the points needed to evaluate the steady-state,
i.e. at late times. Used to optimize memory management and speed of the solver. Default value is an empty list, `saveat=[]`, resulting
in the solution saved at optimal time steps withing the entire time span.

# Examples
```julia-repl
julia> diffeqsolver(x0, p0, tspan, J, bfields, matrix; saveat=((N*4÷5):1:N)*Δt)
```
"""
function diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])

    u0 = [x0[1], x0[2], x0[3], p0[1], p0[2], p0[3], 0, 0, 0, 0, 0, 0]
    Cω2 = matrix.C*transpose(matrix.C)
    bn = t -> matrix.C*[bfields[1](t), bfields[2](t), bfields[3](t)];
    
    function f(du, u, par, t)
        x = @view u[1:3] # @view does not allocate values. No hard copy, just reference.
        p = @view u[4:6]
        v = @view u[7:9]
        w = @view u[10:12]
        du[1:3] = p
        du[4:6] = -(Ω^2)*x + bn(t) + Cω2*v
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
