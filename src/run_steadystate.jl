include("./SpiDy.jl")
using .SpiDy
using NPZ
using ProgressMeter
using Random
using Statistics
using LinearAlgebra

rescaling = false # use if you want to compare to python code
if rescaling
    cfac = 10*(1.76E11)*(1.05E-34)/(1.38E-23)/2
else
    cfac = 1
end

Δt = 0.15
N = 100_000
tspan = (0., N*Δt)
saveat = ((N*4÷5):1:N)*Δt

# Lorentzian(α, ω0, Γ)
J = LorentzianSD(1., 7., 5.); # prm 5
# J = LorentzianSD(100., 7., 5.); # prm 9
matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);
T = 10 .^ LinRange(-3, 3, 24) / cfac

navg = 24 # number of stochastic field realizations to average
p = Progress(length(T));
Sss = zeros(length(T), 3)

for n in 1:length(T)
    noise = ClassicalNoise(T[n]);
    s = zeros(navg, 3)
    Threads.@threads for i in 1:navg
        s0 = [0., 0., -1.] #normalize(rand(3))
        bfields = [bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, matrix; saveat=saveat);
        s[i, :] = mean(sol[2], dims=1)
    end
    Sss[n, :] = mean(s, dims=1)
    next!(p)
end

npzwrite("./notebooks/test.npz", Dict("T" => T*cfac, "S" => Sss))