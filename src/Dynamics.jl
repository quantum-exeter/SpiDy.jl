"""
```Julia
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1], saveat=[])
```

Returns `[sol.t, s, sinterp]`, that is, the vector `sol.t` of time steps at which the solutions are evaluated,
the 3-vector of the solutions `s[1]`, `s[2]`, `s[3]` evaluated at times `sol.t`,
the 3 functions sinterp(t)[i] interpolations of the solutions `s[i]` found in the give time span.

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
    
    function f(du, u, p, t)
        s = @view u[1:3] # @view does not allocate values. No hard copy, just reference.
        v = @view u[4:6]
        w = @view u[7:9]
        du[1:3] = cross(s, Bext + bn(t)/S0 + Cω2*v)
        du[4:6] = w
        du[7:9] = -J.ω0^2*v -J.Γ*w +J.α*s
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