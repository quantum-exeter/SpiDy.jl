using SpiDy
using DifferentialEquations
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
N = 100_000
tspan = (0., N*Δt)
saveat = ((N*4÷5):1:N)*Δt # save diffeqsolver solution at this points
α = 10.
ω0 = 7.
Γ = 5.
J = LorentzianSD(α, ω0, Γ) # coloring the noise
matrix = AnisoCoupling([-sin(π/4) 0. 0. # coupling to the environment
                        0. 0. 0.
                        cos(π/4) 0. 0.]);
T = 10 .^ LinRange(-3, 3, 14) # temperature scale
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
progress = Progress(length(T));
Sss = zeros(length(T), 3*nspin)
for n in eachindex(T)
    noise = QuantumNoise(T[n]);
    s = zeros(navg, 3*nspin)
    Threads.@threads for i in 1:navg
        bfields = [bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, matrix; JH=JH, saveat=saveat, alg=Rodas4());
        s[i, :] = mean(sol, dims=2)
    end
    Sss[n, :] = mean(s, dims=1)
    next!(progress)
end

########################
########################

### Save data NPZ ###
npzwrite("./ss.npz", Dict("T" => T,
                          "S" => Sss))

### Save data CSV ###
dataframe = DataFrame(T = T,
                      Sssx = Sss[:, 1],
                      Sssy = Sss[:, 2],
                      Sssz = Sss[:, 3])
CSV.write("./steadystate.csv", dataframe)

### Plots ###
plot(T, Sss[:, 1], xscale=:log10, xlabel="T", ylabel="S_x")
savefig("./sssx.pdf")

plot(T, Sss[:, 2], xscale=:log10, xlabel="T", ylabel="S_y")
savefig("./sssy.pdf")

plot(T, Sss[:, 3], xscale=:log10, xlabel="T", ylabel="S_z")
savefig("./sssz.pdf")