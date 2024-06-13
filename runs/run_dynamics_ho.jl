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

Δt = 0.01
N = 50_000
tspan = (0., N*Δt)
saveat = (0:1:N)*Δt

J = LorentzianSD(10., 7., 5.); # (α, ω0, Γ)

matrix = IsoCoupling(1.);

noise = ClassicalNoise(10.);

x0 = [1., 0., 0.]
p0 = [0., 0., 0.]

nosc = div(length(x0), 3)

navg = 10

########################
########################

println("Starting...")

progress = Progress(navg);
solx = zeros(navg, length(saveat), 3*nosc)
solp = zeros(navg, length(saveat), 3*nosc)

Threads.@threads for i in 1:navg
    bfields = [bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise),
               bfield(N, Δt, J, noise)];
    sol = diffeqsolver(x0, p0, tspan, J, bfields, matrix; saveat=saveat);
    solx[i, :, :] = transpose(Array(sol)[1:3*nosc, :])
    solp[i, :, :] = transpose(Array(sol)[1+3*nosc:6*nosc, :])
    next!(progress)
end

solavgx = mean(solx, dims=1)[1, :, :];
solavgp = mean(solp, dims=1)[1, :, :];

########################
########################

### Save data NPZ ###
npzwrite("./dynamics.npz", Dict("t" => saveat,
                                "x" => solavgx,
                                "p" => solavgp))

### Save data CSV ###
dataframe = DataFrame(t = saveat,
                      xx = solavgx[:, 1],
                      xy = solavgx[:, 2],
                      xz = solavgx[:, 3],
                      px = solavgp[:, 1],
                      py = solavgp[:, 2],
                      pz = solavgp[:, 3])
CSV.write("./dynamics.csv", dataframe)

### Plots ###
plot(saveat, solavgx[:, 1], xlabel="t", ylabel="x_x")
savefig("./xx.pdf")

plot(saveat, solavgx[:, 2], xlabel="t", ylabel="x_y")
savefig("./xy.pdf")

plot(saveat, solavgx[:, 3], xlabel="t", ylabel="x_z")
savefig("./xz.pdf")

plot(saveat, solavgp[:, 1], xlabel="t", ylabel="p_x")
savefig("./px.pdf")

plot(saveat, solavgp[:, 2], xlabel="t", ylabel="p_y")
savefig("./py.pdf")

plot(saveat, solavgp[:, 3], xlabel="t", ylabel="p_z")
savefig("./pz.pdf")
