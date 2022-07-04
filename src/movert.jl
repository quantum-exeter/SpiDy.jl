include("./SpiDy.jl")
using .SpiDy

using Statistics

using ProgressBars
using NPZ

Δt = 0.15
N = 72_000
tspan = (0., N*Δt)

J = LorentzianSD(1., 7., 5.);

matrix = AnisoCoupling([-sin(π/4) 0. 0.
          0. 0. 0.
          cos(π/4) 0. 0.]);

T = 10 .^ LinRange(-3, 2, 12)
Sss = zeros(length(T), 3)

navg = 12

Threads.@threads for n in ProgressBar(1:length(T))
    noise = ClassicalNoise(T[n]);
    s = zeros(navg, 3)
    for k in 1:navg
        bfields = [bfield(N, Δt, J, noise),
                bfield(N, Δt, J, noise),
                bfield(N, Δt, J, noise)];
        sol = diffeqsolver([0.8, 0., -0.6], tspan, J, bfields, matrix);
        # s[k,:] = mean(sol[3][end-(length(sol[2])÷4):end,:], dims=1)
        s[k,:] = mean([sol[1](t)[n] for t in (3*N÷4)*Δt:Δt:N*Δt, n in 1:3], dims=1)
    end
    Sss[n,:] = mean(s, dims=1)
end

npzwrite("run.npz", Dict("T" => T, "S" => Sss))
