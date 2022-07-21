using SpiDy
using NPZ
using ProgressMeter
using Random
using Statistics
using LinearAlgebra
using Plots

########################
########################

Δt = 0.1
N = 10_000
tspan = (0., N*Δt)
saveat = (0:1:N)*Δt

J = LorentzianSD(10., 7., 5.); # (α, ω0, Γ)

matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);

noise = ClassicalNoise(1.);

s0 = [0., 0., -1.] # normalize(rand(3))

navg = 4

########################
########################

p = Progress(navg);
sols = zeros(navg, length(saveat), 3)
Threads.@threads for i in 1:navg
    bfields = [bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise)];
    sol = diffeqsolver(s0, tspan, J, bfields, matrix; saveat=saveat);
    sols[i, :, :] = sol[2]
    next!(p)
end
solavg = mean(sols, dims=1)[1, :, :];

########################
########################

npzwrite("./dynamics.npz", Dict("t" => saveat, "S" => solavg))

plot(saveat, solavg[:, 1], xlabel="t", ylabel="S_x")
savefig("./sx.pdf")

plot(saveat, solavg[:, 2], xlabel="t", ylabel="S_y")
savefig("./sy.pdf")

plot(saveat, solavg[:, 3], xlabel="t", ylabel="S_z")
savefig("./sz.pdf")