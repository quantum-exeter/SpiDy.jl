include("./SpiDy.jl")
using .SpiDy

using Statistics
using ProgressBars
using NPZ
using Random
using LinearAlgebra

Δt = 0.15
N = 72_000
tspan = (0., N*Δt)

# Lorentzian(α, ω0, Γ)
J = LorentzianSD(1., 7., 5.); # prm 5

matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);

T = 10 .^ LinRange(-2, 2, 12)
Sss = zeros(length(T), 3)

navg = 5
Threads.@threads for n in ProgressBar(1:length(T))
    noise = ClassicalNoise(T[n]);
    s = zeros(navg, 3)
    for i in 1:navg
        s0 = normalize(rand(3)) # [0.8, 0., -0.6]
        bfields = [bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, matrix);
        s[i, :] = mean([sol[3](t)[n] for t in (N*3÷4)*Δt:Δt:N*Δt, n in 1:3], dims=1)
    end
    Sss[n, :] = mean(s, dims=1)
end

npzwrite("../notebooks/data.npz", Dict("T" => T, "S" => Sss))