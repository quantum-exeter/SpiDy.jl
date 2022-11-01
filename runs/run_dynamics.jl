using SpiDy
using NPZ
using DataFrames
using CSV
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
α = 10.
ω0 = 7.
Γ = 5.
J = LorentzianSD(α, ω0, Γ) # coloring the noise
matrix = AnisoCoupling([-sin(π/4) 0. 0. # coupling to the environment
                        0. 0. 0.
                        cos(π/4) 0. 0.]);
T = 1.
noise = ClassicalNoise(T);
navg = 6 # number of stochastic realizations
nspin = 4 # number of spins
s0 = zeros(3*nspin)
for i in 1:nspin
    ϵ = 0.1
    s0tmp = [ϵ*rand(), ϵ*rand(), -1]
    s0[1+(i-1)*3:3+(i-1)*3] = s0tmp./norm(s0tmp)
end
J0 = 1.
JH = Nchain(nspin, J0)

########################
########################

println("Starting...")
progress = Progress(navg);
sols = zeros(navg, length(saveat), 3*nspin)
Threads.@threads for i in 1:navg
    bfields = [bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise)];
    sol = diffeqsolver(s0, tspan, J, bfields, matrix; JH=JH, saveat=saveat);
    sols[i, :, :] = transpose(sol[:, :])
    next!(progress)
end
solavg = mean(sols, dims=1)[1, :, :];

########################
########################

### Save data NPZ ###
npzwrite("./dynamics.npz", Dict("t" => saveat,
                                "S" => solavg))

### Save data CSV ###
dataframe = DataFrame(t = saveat,
                      Sx = solavg[:, 1],
                      Sy = solavg[:, 2],
                      Sz = solavg[:, 3])
CSV.write("./dynamics.csv", dataframe)

### Plots ###
plot(saveat, solavg[:, 1], xlabel="t", ylabel="S_x")
savefig("./sx.pdf")

plot(saveat, solavg[:, 2], xlabel="t", ylabel="S_y")
savefig("./sy.pdf")

plot(saveat, solavg[:, 3], xlabel="t", ylabel="S_z")
savefig("./sz.pdf")
