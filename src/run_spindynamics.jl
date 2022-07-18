include("./SpiDy.jl")
using .SpiDy
using NPZ
using ProgressMeter
using Random
using Statistics
using LinearAlgebra

Δt = 0.15 *2π
N = 1_000
tspan = (0., N*Δt)
saveat = ((N*4÷5):1:N)*Δt

# Lorentzian(α, ω0, Γ)
# J = LorentzianSD(1., 7., 5.); # prm 5
J = LorentzianSD(100., 7., 5.); # prm 9

matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);

T = 10 .^ LinRange(-3, 3, 12)
Sss = zeros(length(T), 3)

noise = ClassicalNoise(1.);
s0 = [0.8, 0., 0.6] #normalize(rand(3))
bfields = [bfield(N, Δt, J, noise),
           bfield(N, Δt, J, noise),
           bfield(N, Δt, J, noise)];
           
sol = diffeqsolver(s0, tspan, J, bfields, matrix);

npzwrite("./notebooks/data_dynamics_prm9_b(t)=0.npz", Dict("T" => sol[1], "S" => sol[2]))