using SpiDy
using NPZ
using ProgressMeter
using Random
using Statistics
using LinearAlgebra
using Plots

########################
########################

rescaling = false # use if you want to compare to python code
if rescaling
    cfac = 10*(1.76E11)*(1.05E-34)/(1.38E-23)/2
else
    cfac = 1
end

Δt = 0.1
N = 100_000
tspan = (0., N*Δt)
saveat = ((N*4÷5):1:N)*Δt

J = LorentzianSD(10., 7., 5.); # (α, ω0, Γ)

matrix = AnisoCoupling([-sin(π/4) 0. 0.
                        0. 0. 0.
                        cos(π/4) 0. 0.]);

T = 10 .^ LinRange(-3, 3, 12) / cfac

navg = 6 # number of stochastic field realizations to average

########################
########################

println("Starting...")

p = Progress(length(T));
Sss = zeros(length(T), 3)

for n in 1:length(T)
    noise = ClassicalNoise(T[n]);
    s = zeros(navg, 3)
    Threads.@threads for i in 1:navg
        s0 = [0., 0., -1.] # normalize(rand(3))
        bfields = [bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, matrix; saveat=saveat);
        s[i, :] = mean(sol[2], dims=1)
    end
    Sss[n, :] = mean(s, dims=1)
    next!(p)
end

########################
########################

npzwrite("./steadystate.npz", Dict("T" => T*cfac, "S" => Sss))

plot(T, Sss[:, 1], xscale=:log10, xlabel="T", ylabel="S_x")
savefig("./sssx.pdf")

plot(T, Sss[:, 2], xscale=:log10, xlabel="T", ylabel="S_y")
savefig("./sssy.pdf")

plot(T, Sss[:, 3], xscale=:log10, xlabel="T", ylabel="S_z")
savefig("./sssz.pdf")