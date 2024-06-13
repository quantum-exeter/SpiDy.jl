# We show an example code to generate a run of SpiDy.jl for a single spin interacting with an environment.
# Given the stochastic nature of the problem solved, we will be dealing with different dynamical trajectories.
# In the following code, we show the parameters needed to build a subset of these trajectories, solutions
# to the stochastic differential equations of Eq.(1). Then, we plot a single one of them as an example.
# The entire code is commented throughout for a better understanding of the single elements of the run.
# We show the results of the dynamics averaged over a larger set (10000) of trajectories in \autoref{fig:dynamics}.
# Note that both the average and the standard deviation of the set of trajectories are not evaluated with the
# following code but are nonetheless represented in the figure for clarity. Further examples are available in the code repository.


### importing SpiDy ###
using SpiDy

### setting the parameters ###
ωL = 1 # Larmor frequency (reference time scale)
Δt = 0.1 / ωL # time step for the dynamics evaluation
tend = 150 / ωL # final time of the dynamics
N = round(Int, tend/Δt) # number of total steps
tspan = (0, N*Δt) # tuple of initial and final time
saveat = (0:1:N)*Δt # vector of times at which the solution is saved
α = 10 * ωL # Lorentzian coupling amplitude
ω0 = 7 * ωL # Lorentzian resonant frequency
Γ = 5 * ωL # Lorentzian width
Jsd = LorentzianSD(α, ω0, Γ) # Lorentzian spectral density
Cw = IsoCoupling(1) # isotropic coupling tensor
# the resulting coupling tensor is equivalent to the following
# Cw = AnisoCoupling([1 0 0
#                     0 1 0
#                     0 0 1]);
T = 0.8 * ωL # temperature at which the dynamics takes place (where ħ=1, kB=1)
noise = ClassicalNoise(T) # noise profile for the stochastic field
s0 = [1.0; 0.0; 0.0] # initial conditions of the spin vector for the dynamics
ntraj = 10 # number of trajectories (stochastic realizations)

### running the dynamics ###
sols = zeros(ntraj, 3, length(saveat)) # solution matrix
for i in 1:ntraj # iterations through the number of trajectories
    # we use the Lorentzian spectral density Jsd to generate the stochastic
    # field. This ensures the field obeys the FDR as noted in the main text
    local bfields = [bfield(N, Δt, Jsd, noise),
                     bfield(N, Δt, Jsd, noise), # vector of independent
                     bfield(N, Δt, Jsd, noise)] # stochastic fields
    # diffeqsolver (below) solves the system for the single trajectory
    local sol = diffeqsolver(s0, tspan, Jsd, bfields, Cw; saveat=saveat)
    sols[i,:,:] = Array(sol) # store the trajectory into the matrix of solutions
end

### example plot ###
# use Plots.jl pkg to plot a single trajectory of the dynamics over time
using Plots
plot(xlabel="time", ylabel="spin components")
# sols[i,j,k] with i: trajectory index, j: spin component, k: solution at
# the k-th time point
plot!(saveat, sols[1,1,:], label="x-component")
plot!(saveat, sols[1,2,:], label="y-component")
plot!(saveat, sols[1,3,:], label="z-component")
savefig("example_run.pdf")
