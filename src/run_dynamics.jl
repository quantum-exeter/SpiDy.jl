include("./SpiDy.jl")
using .SpiDy
using NPZ
using ProgressMeter
using Random
using Statistics
using LinearAlgebra

################

Δt = 0.00015
N = 1000_000
tspan = (0., N*Δt)
saveat = (0:1:N)*Δt

# Lorentzian(α, ω0, Γ)
J = LorentzianSD(1., 7., 5.); # prm 5
# J = LorentzianSD(100., 7., 5.); # prm 9

matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);
noise = ClassicalNoise(0.1);
s0 = [0., 0., -1.] #normalize(rand(3))

navg = 10

###############

p = Progress(navg);
sols = zeros(navg, length(saveat), 3)
for i in 1:navg
    bfields = [bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise)];
    sol = diffeqsolver(s0, tspan, J, bfields, matrix; saveat=saveat);
    sols[i, :, :] = sol[2]
    next!(p)
end
solavg = mean(sols, dims=1)[1, :, :];

npzwrite("./notebooks/test3.npz", Dict("T" => saveat, "S" => solavg))