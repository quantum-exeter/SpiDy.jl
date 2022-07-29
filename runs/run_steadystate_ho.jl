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
N = 50_000
tspan = (0., N*Δt)
saveat = ((N*4÷5):1:N)*Δt

J = LorentzianSD(10., 7., 5.); # (α, ω0, Γ)

matrix = IsoCoupling(1.);

T = 10 .^ LinRange(-3, 3, 12)

navg = 10 # number of stochastic field realizations to average

########################
########################

println("Starting...")

progress = Progress(length(T));
xss = zeros(length(T), 3)
pss = zeros(length(T), 3)

for n in eachindex(T)
    noise = ClassicalNoise(T[n]);
    x = zeros(navg, 3)
    p = zeros(navg, 3)
    Threads.@threads for i in 1:navg
        x0 = [1., 0., 0.]
        p0 = [0., 0., 0.]
        bfields = [bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise),
                   bfield(N, Δt, J, noise)];
        sol = diffeqsolver(x0, p0, tspan, J, bfields, matrix; saveat=saveat);
        x[i, :] = mean(sol[2].^2, dims=1)
        p[i, :] = mean(sol[3].^2, dims=1)
    end
    xss[n, :] = mean(x, dims=1)
    pss[n, :] = mean(p, dims=1)
    next!(progress)
end

########################
########################

### Save data NPZ ###
# npzwrite("./steadystate.npz", Dict("T" => T,
#                                    "x" => xss,
#                                    "p" => pss))

### Save data CSV ###
# dataframe = DataFrame(T = T,
#                       xssx = xss[:, 1],
#                       xssy = xss[:, 2],
#                       xssz = xss[:, 3],
#                       pssx = pss[:, 1],
#                       pssy = pss[:, 2],
#                       pssz = pss[:, 3])
# CSV.write("./steadystate.csv", dataframe)

### Plots ###
plot(T, xss[:, 1], xscale=:log10, xlabel="T", ylabel="x_x")
savefig("./xssx.pdf")

plot(T, xss[:, 2], xscale=:log10, xlabel="T", ylabel="x_y")
savefig("./xssy.pdf")

plot(T, xss[:, 3], xscale=:log10, xlabel="T", ylabel="x_z")
savefig("./xssz.pdf")

plot(T, pss[:, 1], xscale=:log10, xlabel="T", ylabel="p_x")
savefig("./pssx.pdf")

plot(T, pss[:, 2], xscale=:log10, xlabel="T", ylabel="p_y")
savefig("./pssy.pdf")

plot(T, pss[:, 3], xscale=:log10, xlabel="T", ylabel="p_z")
savefig("./pssz.pdf")