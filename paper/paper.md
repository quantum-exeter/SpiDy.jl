---
title: 'SpiDy.jl: open-source Julia package for the study of non-Markovian stochastic dynamics'
tags:
  - Julia
  - spin-boson model
  - harmonic oscillator
  - stochastic field
  - non-Markovian memory
  - colored noise
  - anisotropic coupling
  - magnetism
authors:
  - name: Stefano Scali
    orcid: 0000-0002-8133-1551
    corresponding: true
    affiliation: 1
  - name: Simon Horsley
    orcid: 0000-0001-7242-7941
    affiliation: 1 
  - name: Janet Anders
    orcid: 0000-0002-9791-0363
    affiliation: "1, 2"
  - name: Federico Cerisola
    orcid: 0000-0003-2961-739X
    affiliation: "1, 3"
affiliations:
 - name: Department of Physics and Astronomy, University of Exeter, Exeter EX4 4QL, United Kingdom
   index: 1
 - name: Institute of Physics and Astronomy, University of Potsdam, 14476 Potsdam, Germany
   index: 2
 - name: Department of Engineering Science, University of Oxford, Oxford OX1 3PJ, United Kingdom.
   index: 3
date: 4 October 2023
bibliography: paper.bib
---

# Summary

SpiDy.jl solves the non-Markovian stochastic dynamics of interacting classical spin vectors and harmonic oscillator networks in contact with a dissipative environment.
The methods implemented allow the user to include arbitrary memory effects and colored quantum noise spectra.
In this way, SpiDy.jl provides key tools for the simulation of classical and quantum open systems including non-Markovian effects and arbitrarily strong coupling to
the environment. Among the wide range of applications, some examples range from atomistic spin dynamics to ultrafast magnetism and the study of anisotropic materials.
We provide the user with Julia notebooks to guide them through the various mathematical methods and help them quickly set up complex simulations.

# Statement of need

The problem of simulating the dynamics of interacting rotating bodies and harmonic oscillator networks in the presence of a dissipative environment can find a vast range of applications in the modeling of physical systems. This task is rendered particularly challenging when one desires to capture the non-Markovian effects that arise in the dynamics due to strong coupling with the environment. SpiDy.jl is a library that allows the user to efficiently simulate these systems to obtain both detailed dynamics and steady-state properties.

A relevant example of the applicability of SpiDy.jl is the modeling of spins at low temperatures and at short timescales, which is a fundamental task to address many open questions in the field of magnetism and magnetic material modeling [@halilov98].
State-of-the-art tools such as those developed for atomistic spin dynamics simulations are based on solving the Landau–Lifshitz–Gilbert (LLG) equation [@evans2014].
Despite their massive success, these tools run into shortcomings in accurately modeling systems at low temperatures and for short timescales where environment memory effects have been observed [@ciornei2011; @neeraj2020].
Recent work has focused on developing a comprehensive quantum-thermodynamically consistent framework suitable to model the dynamics of spins in magnetic materials while addressing these shortcomings [@Anders_2022]. This framework includes strong coupling effects to the environment such as non-Markovian memory, colored noise, and quantum-like fluctuations.
At its core, SpiDy.jl implements the theoretical framework introduced in [@Anders_2022], allowing for the study of environment memory effects and anisotropic system-environment coupling.
SpiDy.jl can be readily adopted for atomistic spin dynamics simulations [@evans2014; @Barker_2019], ultrafast magnetism [@Beaurepaire_1996], and ferromagnetic and semiconductive systems exhibiting anisotropic damping [@Chen_2018].
A further set of applications stems from the extension of SpiDy.jl to handle the non-Markovian stochastic dynamics of harmonic oscillators. This model will be of interest in the field of quantum thermodynamics where harmonic oscillators play a key role in modeling open quantum systems.
The package is written in pure Julia to take advantage of the language performance.

The software package has seen a wide range of applications to date. Firstly, the convenience of three independent environments in SpiDy.jl finds application in the microscopic modeling of spins affected by noise due to vibrations of the material lattice [@Anders_2022]. SpiDy.jl also found application in the demonstration of the quantum-to-classical correspondence at all coupling strengths between a spin and an external environment [@cerisola2022quantumclassical]. Here, the temperature dependence of the spin steady-state magnetization obtained with SpiDy.jl is successfully compared with the classical mean force state of the system. In Ref. [@hartmann2023anisotropic], the authors take advantage of the customizable coupling tensor in SpiDy.jl to explore the anisotropic effects of the environment on the system. In Ref. [@berritta2023accounting], SpiDy.jl is used as a sub-routine to build quantum-improved atomistic spin dynamics simulations. In the paper, the authors take advantage of the customizable power spectrum to implement ad-hoc simulations matching known experimental results. Lastly, with an eye to the harmonic oscillator side, SpiDy.jl is used to match the quantum harmonic oscillator dynamics with its stochastic counterpart [@glatthard2023harmonic]. Here, the authors exploit the recent implementation of harmonic oscillator dynamics.

# Overview

To model a system of interacting classical spin vectors, SpiDy.jl solves the generalized stochastic LLG equation [@Anders_2022]
$$
\frac{\mathrm{d} \mathbf{S}_n(t)}{\mathrm{d} t} = \frac{1}{2} \mathbf{S}_n(t) \times
\bigg[\sum_{m\neq n}J_{n,m}\mathbf{S}_m(t) + \mathbf{B} + \mathbf{b}_n(t) + \int_{t_0}^t\mathrm{d}t^{\prime} K_n(t-t^{\prime}) \mathbf{S}_n(t^{\prime})\bigg], \quad (1)
$$
where $\mathbf{S}_n(t)$ represents the $n$-th spin vector, the interaction matrix $J_{n,m}$ sets the interaction strength between the $n$-th and $m$-th spins, $\mathbf{B}$ is the external field, which determines the natural precession direction and frequency of the spins in the absence of interaction, and $b_n(t)$ is the time-dependent stochastic field induced by the environment. Finally, the last integral term in Eq.(1) gives the spin dissipation due to the environment, including non-Markovian effects accounted for by the memory kernel matrix $K_n(t)$.
Here, we allow each spin to interact with three independent sources of noise so that in general $K_n(t) = C_n k_n(t)$, where $k_n(t)$ is a time-dependent function and $C_n$ is a $3\times3$ matrix that determines how each of the $n$-th spin components couples to each of the three noise sources.
To efficiently simulate the non-Markovian effects, we follow the methods explained in [@Anders_2022] and restrict ourselves to the case where the memory kernel $k(t)$ comes from a Lorentzian spectral density of the bath $\mathcal{J}(\omega) = \alpha\Gamma/((\omega_0^2 - \omega^2)^2 + \Gamma^2\omega^2)$ with peak frequency $\omega_0$, peak width $\Gamma$ and amplitude $\alpha$, so that $k(t) = \Theta(t)\,\alpha\,e^{-\Gamma t/2}\sin(\omega_1 t)/\omega_1$, where $\omega_1^2 = \omega_0^2 - \Gamma^2/4$.
In the code, the stochastic noise $b_n(t)$ is generated so that it satisfies the fluctuation-dissipation relation (FDR) (see [@Anders_2022]). That is, the power spectral density of the stochastic noise satisfies $P(\omega, T) = \mathcal{J}(\omega) \eta(T)$ where $\mathcal{J}(\omega)$ is the Lorentzian spectral density and $\eta(T)$ defines the temperature dependence on the bath. Here, the user can choose, among others, a classical or quantum-like temperature dependence, namely $\eta_\mathrm{cl}(T) = k_\mathrm{B}T/2\hbar\omega$ and $\eta_\mathrm{qu}(T) = \mathrm{coth}(2\hbar\omega/k_\mathrm{B}T)$ respectively.

In addition, SpiDy.jl also allows one to study the stochastic dynamics of coupled harmonic oscillator networks. In the same way, as for the spin case, the harmonic oscillators can be coupled together with a user-defined system-system interaction. The equations of motion solved in this case are
$$
\frac{\mathrm{d}^2 \mathbf{X}_n(t)}{\mathrm{d} t^2} =
\sum_{m\neq n}J_{n,m}\mathbf{X}_m(t) - \Omega^2 \mathbf{X}_n(t) + \mathbf{b}_n(t) + \int_{t_0}^t\mathrm{d}t^{\prime} K_n(t-t^{\prime}) \mathbf{X}_n(t^{\prime}), \quad (2)
$$
where $\mathbf{X}_n(t)$ represents the position vector of the $n$-th harmonic oscillator, the interaction matrix $J_{n,m}$ sets the interaction strength between the $n$-th and $m$-th harmonic oscillators, and $\Omega$ is the bare frequency of the harmonic oscillators (we consider identical oscillators). All other terms have the same role as in the spin case (see Eq.(1)).

In conclusion, SpiDy.jl implements the stochastic dynamics of coupled integro-differential equations to model systems of interacting spins or harmonic oscillator networks subject to environment noise.
Among others, some of the key features of the package include:

- Coloured stochastic noise that satisfies the FDR and accounts for both classical and quantum bath statistics.
- Simulation of non-Markovian system dynamics due to the memory kernel $K_n(t)$.
- Custom system-environment coupling tensors, allowing for isotropic or anisotropic couplings. Both amplitudes and geometry of the coupling can be specified.
- Choice between local environments, i.e. distinct baths acting on the single sub-system, or a single common environment.

In the next section, we show a minimal working example to run the spin dynamics where we list all the required input parameters.

# Example

![**Single-spin dynamics.** Dynamics of the $x$, $y$, and 
 $z$ spin components. The components are normalized against the total spin length $S_0$ and time axes are expressed in units of the Larmor frequency $\omega_\mathrm{L}$ ($\omega_\mathrm{L} = |\mathbf{B}|$ in Eq.(1)). We show an example set of 5 stochastic trajectories of the spin dynamics (colored semi-transparent lines) together with their stochastic average (gray solid line). Note that, while we show only 5 trajectories for clarity, the average dynamics is obtained from 10000 trajectories. We also represent the range of one standard deviation from the average dynamics (gray-shaded area). In the inset, we show the convergence of the same dynamics towards the steady state at longer times. This example is obtained using the Lorentzian parameters "set 1" found in Ref. [@Anders_2022]. The code used to generate the stochastic trajectories is shown in the text. \label{fig:spin_dynamics}](spin_dynamics.pdf){ width=100% }

Now, we show an example code to generate a run of SpiDy.jl for a single spin interacting with an environment. Given the stochastic nature of the problem solved, we will be dealing with different dynamical trajectories. In the following code, we show the parameters needed to build a subset of these trajectories, solutions to the stochastic differential equations of Eq.(1). Then, we plot a single one of them as an example. The entire code is commented throughout for a better understanding of the single elements of the run. We show the results of the dynamics averaged over a larger set (10000) of trajectories in \autoref{fig:spin_dynamics}. Note that both the average and the standard deviation of the set of trajectories are not evaluated with the following code but are nonetheless represented in the figure for clarity. Further examples are available in the code repository.

```Julia
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
    sols[i,:,:] = sol[:,:] # store the trajectory into the matrix of solutions
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
```

# Acknowledgements

SS and FC thank Marco Berritta and Charlie Hogg for insightful suggestions for the implementation of the spin-spin coupling and the harmonic oscillator dynamics.
SARH and JA thank the Royal Society for their support. SS is supported by a DTP grant from EPSRC (EP/R513210/1). SARH acknowledges the Royal Society and TATA for financial support through the Grant URFR 211033. JA and FC acknowledge funding from EPSRC (EP/R045577/1). FC gratefully acknowledges funding from the Foundational Questions Institute Fund (FQXi–IAF19-01).

# References
