var documenterSearchIndex = {"docs":
[{"location":"couplingfunctions/#Spin-spin-interactions","page":"Spin-spin interactions","title":"Spin-spin interactions","text":"","category":"section"},{"location":"couplingfunctions/","page":"Spin-spin interactions","title":"Spin-spin interactions","text":"Nchain\nNNlattice","category":"page"},{"location":"couplingfunctions/#SpiDy.Nchain","page":"Spin-spin interactions","title":"SpiDy.Nchain","text":"Nchain(N, J0; boundary=nothing)\n\nCreate the spin-spin coupling matrix of a 1D chain of N sites with nearest-neighbour interactions of strength set by J0.\n\nArguments\n\nN: The number of sites in the chain.\nJ0: The nearest-neighbour coupling strength.\nboundary=nothing: (Optional) Specifies the boundary condition of the chain. Default is nothing which corresponds to a chain open at the edges. Use :periodic for periodic boundary conditions.\n\nReturns\n\nAn N×N array representing the coupling matrix.\n\nExamples\n\njulia> J = Nchain(4, 2.0)\n4×4 Matrix{Float64}:\n 0.0  2.0  0.0  0.0\n 2.0  0.0  2.0  0.0\n 0.0  2.0  0.0  2.0\n 0.0  0.0  2.0  0.0\n\njulia> J_periodic = Nchain(4, 2.0; boundary=:periodic)\n4×4 Matrix{Float64}:\n 0.0  2.0  0.0  2.0\n 2.0  0.0  2.0  0.0\n 0.0  2.0  0.0  2.0\n 2.0  0.0  2.0  0.0\n\n\n\n\n\n","category":"function"},{"location":"couplingfunctions/#SpiDy.NNlattice","page":"Spin-spin interactions","title":"SpiDy.NNlattice","text":"NNlattice(N, Jh, Jv; boundary=nothing)\n\nCreate the spin-spin coupling matrix of a NxN 2D lattice with nearest-neighbour interactions and specified horizontal and vertical coupling strengths, Jh and Jv respectively.\n\nArguments\n\nN: The size of the lattice (N x N).\nJh: The horizontal nearest-neighbour coupling strength.\nJv: The vertical nearest-neighbour coupling strength.\nboundary=nothing: (Optional) Specifies the boundary condition of the lattice. Default is nothing which corresponds to open edges. Use :periodic for periodic boundary condition.\n\nReturns\n\nAn N^2×N^2 array representing the coupling matrix of the lattice.\n\nExamples\n\njulia> J = NNlattice(3, 2.0, 1.0)\n9×9 LinearAlgebra.Symmetric{Float64, Matrix{Float64}}:\n 0.0  2.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0\n 2.0  0.0  2.0  0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  2.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0\n 1.0  0.0  0.0  0.0  2.0  0.0  1.0  0.0  0.0\n 0.0  1.0  0.0  2.0  0.0  2.0  0.0  1.0  0.0\n 0.0  0.0  1.0  0.0  2.0  0.0  0.0  0.0  1.0\n 0.0  0.0  0.0  1.0  0.0  0.0  0.0  2.0  0.0\n 0.0  0.0  0.0  0.0  1.0  0.0  2.0  0.0  2.0\n 0.0  0.0  0.0  0.0  0.0  1.0  0.0  2.0  0.0\n\njulia> J_periodic = NNlattice(2, 2.0, 1.0; boundary=:periodic)\n9×9 LinearAlgebra.Symmetric{Float64, Matrix{Float64}}:\n 0.0  2.0  2.0  1.0  0.0  0.0  1.0  0.0  0.0\n 2.0  0.0  2.0  0.0  1.0  0.0  0.0  1.0  0.0\n 2.0  2.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0\n 1.0  0.0  0.0  0.0  2.0  2.0  1.0  0.0  0.0\n 0.0  1.0  0.0  2.0  0.0  2.0  0.0  1.0  0.0\n 0.0  0.0  1.0  2.0  2.0  0.0  0.0  0.0  1.0\n 1.0  0.0  0.0  1.0  0.0  0.0  0.0  2.0  2.0\n 0.0  1.0  0.0  0.0  1.0  0.0  2.0  0.0  2.0\n 0.0  0.0  1.0  0.0  0.0  1.0  2.0  2.0  0.0\n\n\n\n\n\n","category":"function"},{"location":"noise/#Noise","page":"Noise","title":"Noise","text":"","category":"section"},{"location":"noise/","page":"Noise","title":"Noise","text":"Noise\nClassicalNoise\nQuantumNoise\nNoZeroQuantumNoise\nspectrum\npsd","category":"page"},{"location":"noise/#SpiDy.Noise","page":"Noise","title":"SpiDy.Noise","text":"abstract type Noise\n\nAn abstract type used to represent different kinds of environment noise.\n\nAny user-defined subtype of Noise should implement a method spectrum(::NoiseSubType, ω) which returns the spectrum of the noise at the given frequency ω.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.ClassicalNoise","page":"Noise","title":"SpiDy.ClassicalNoise","text":"struct ClassicalNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent classical noise, with spectrum\n\nmathcalS(omega) = frac2k_mathrmBThbaromega\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.QuantumNoise","page":"Noise","title":"SpiDy.QuantumNoise","text":"struct QuantumNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent quantum noise, with spectrum\n\nmathcalS(omega) = cothleft(frachbaromega2k_mathrmBTright)\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.NoZeroQuantumNoise","page":"Noise","title":"SpiDy.NoZeroQuantumNoise","text":"struct NoZeroQuantumNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent quantum noise with zero-point-fluctuations removed, that is with spectrum\n\nmathcalS(omega) = cothleft(frachbaromega2k_mathrmBTright) - 1\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.spectrum","page":"Noise","title":"SpiDy.spectrum","text":"spectrum(::Noise)\n\nReturn the spectrum of the given noise as a function of frequency.\n\nArguments\n\nn::Noise: The noise type.\n\nReturns\n\nA function that takes a frequency ω and returns the corresponding spectrum value.\n\n\n\n\n\n","category":"function"},{"location":"noise/#SpiDy.psd","page":"Noise","title":"SpiDy.psd","text":"function psd(J::AbstractSD, noise::Noise)\n\nCalculate the Power Spectral Density (PSD) for a given environment spectral density J and noise model noise.\n\nArguments\n\nJ::AbstractSD: The environment spectral density.\nnoise::Noise: The noise model for the environment.\n\nNote: The LorentzianSD type is provided by the SpectralDensities.jl package.\n\nReturns\n\nA function psd(ω) that calculates the PSD at a given frequency ω.\n\n\n\n\n\n","category":"function"},{"location":"units/#Units-and-choice-of-input-parameters","page":"Units and parameter choice","title":"Units and choice of input parameters","text":"","category":"section"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"In this document we describe the choice of units implemented in SpiDy.jl (which follows the conventions utilised in  NJP 24 033020 (2022)) and how to appropriately choose the parameters that the library takes as input.","category":"page"},{"location":"units/#Equations-of-motion","page":"Units and parameter choice","title":"Equations of motion","text":"","category":"section"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"SpiDy.jl is designed to implement the equations of spins in presence of a bath of harmonic oscillators with Lorentzian spectral density, as derived in NJP 24 033020 (2022), which for a single spin mathbfS in presence of an external field mathbfB_mathrmext read","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"fracmathrmdmathbfSmathrmdt =\n    gammamathbfStimesleft(mathbfS_mathrmext + mathbfb + mathbfVright) \nfracmathrmdmathbfVmathrmdt = mathbfW \nfracmathrmdmathbfWmathrmdt = gamma A mathbfS - omega_0^2mathbfV - GammamathbfW","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"where A, omega_0, and Gamma parametrise the Lorentzian spectral density as","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"J(omega) = fracAGammapi fracomega(omega_0^2 - omega^2)^2 + omega^2Gamma^2","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"and the stochastic mathbfb is given by","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"mathbfb(t) = int_-infty^+inftymathrmdt F(t-t) xi(t)","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"with xi being white noise and","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"F(tau) = frac12piint_-infty^+inftymathrmdomega\n    e^-iomegatau sqrtP(omega)","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"Here, P(omega) is the power spectral density of the environment and it is given in terms of the Lorentzian spectral density J(omega) and the environment thermal noise N(omega) by P(omega) = hbarpi J(omega) N(omega). The noise N(omega) can be classical, quantum, or quantum with no zero point fluctuations. For example, for the quantum case we have","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"N_mathrmqu(omega) = cothleft(frachbaromega2k_mathrmBTright)","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"Note that in these equations above, all quantities have standard units.","category":"page"},{"location":"units/#Units-rescaling","page":"Units and parameter choice","title":"Units rescaling","text":"","category":"section"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"SpiDy.jl implements these equations in a unit-free way, following the conventions of NJP 24 033020 (2022). In summary, suppose a spin of length hbar S_0 is in presence of a magnetic field of magnitude B_0. We define the Larmor frequency","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"omega_mathrmL = gamma_e B_0","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"with gamma_e the electron gyromagnetic ratio.","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"SpiDy.jl then takes as input parameters for the simulation the following unit free parameters:","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"S_0: the spin length (in units of hbar).\nbarB_mathrmext: external magnetic field.\nbart_mathrmend: final time of the evolution.\nmathrmdbart: time differential.\nbaromega_0: peak frequency of the Lorentzian spectral density.\nbarGamma: with of the Lorentzian spectral density.\nbaralpha: amplitude of the Lorentzian density.\nbarT: the environment temperature.","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"These quantities are related to the unitful units in the previous section by the following conversions:","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"beginalign\nmathbfB_mathrmext = B_0  barB_mathrmext \nt_mathrmend = omega_mathrmL^-1  bart_mathrmend \nmathrmdt = omega_mathrmL^-1  mathrmdbart \nomega_0 = omega_mathrmL  baromega_0 \nGamma = omega_mathrmL  barGamma \nomega_0 = omega_mathrmL  baromega_0 \nA = fracB_0^2omega_mathrmLhbar S_0  baralpha \nT = frachbaromega_mathrmLk_mathrmB  barT\nendalign","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"With these definitions, the unit-free Gilbert damping is given by","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"eta = fracbaralphabarGammabaromega_0^4","category":"page"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"warning: Units convetion in old versions of SpiDy.jl\nNote that the unit conversion just explained, in line with NJP 24 033020 (2022), is correct for SpiDy.jl versions 1.2.0 and up.Older version of SpiDy.jl used a different convention, which unfortunately ment that a change of spin length S_0 implied a change in temperature and time-scale, therefore requiring redefinig the paramerts of the Loretnzian, evolution time, etc.The change in units convetion means that versions of SpiDy.jl before and after 1.2.0 will produce different results. This breaking change was made to bring the code into consistency with the article NJP 24 033020 (2022) and to make the mapping of parameters used in SpiDy.jl to real units much easier and more straightforward (see also next section).Finally, it is worth noting that the behaviour of SpiDy.jl in older versions can be exactly recovered by passing the option S_0 = 1 to the integrator, since in that case both units conventions agree.","category":"page"},{"location":"units/#Extracting-Lorentzian-units-from-experimental-data","page":"Units and parameter choice","title":"Extracting Lorentzian units from experimental data","text":"","category":"section"},{"location":"units/","page":"Units and parameter choice","title":"Units and parameter choice","text":"To be written.","category":"page"},{"location":"stochasticfield/#Stochastic-field","page":"Stochastic field","title":"Stochastic field","text":"","category":"section"},{"location":"stochasticfield/","page":"Stochastic field","title":"Stochastic field","text":"bfield","category":"page"},{"location":"stochasticfield/#SpiDy.bfield","page":"Stochastic field","title":"SpiDy.bfield","text":"bfield(N, Δt, J::AbstractSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true)\n\nGenerate a stochastic field realisation time series (of length N and spacing Δt) based on the given noise model noise and spectral density J.\n\nArguments\n\nN: The number of time steps.\nΔt: The time step size.\nJ::AbstractSD: The environment spectral density.\nnoise::Noise: The noise model for the environment.\ndistro=Normal(0., 1/sqrt(Δt)): (Optional) The distribution of noise samples. Default is a normal distribution with mean 0 and standard deviation 1/sqrt(Δt).\ninterpolation=true: (Optional) Specifies whether to use linear interpolation for the stochastic field time series. Default is true.\n\nNote: The AbstractSD type is provided by the SpectralDensities.jl package.\n\nReturns\n\nA time series of the stochastic field values.\n\n\n\n\n\n","category":"function"},{"location":"couplingtensor/#Environment-coupling-tensor","page":"Spin-environment coupling tensor","title":"Environment coupling tensor","text":"","category":"section"},{"location":"couplingtensor/","page":"Spin-environment coupling tensor","title":"Spin-environment coupling tensor","text":"Coupling\nAnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}}\nIsoCoupling{TT<:Real}","category":"page"},{"location":"couplingtensor/#SpiDy.Coupling","page":"Spin-environment coupling tensor","title":"SpiDy.Coupling","text":"abstract type Coupling\n\nAbstract type used to represent different forms of the environment coupling tensor.\n\n\n\n\n\n","category":"type"},{"location":"couplingtensor/#SpiDy.AnisoCoupling","page":"Spin-environment coupling tensor","title":"SpiDy.AnisoCoupling","text":"struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling\n\nA subtype of Coupling used to represent generic anisotropic coupling between the spin and environment.\n\nFields\n\nC::TT: The coupling matrix, which must be a 3x3 matrix (of type Matrix{<:Real}).\n\nExamples\n\nFor a system with environment coupling where the y coupling is twice as large as the x, and the z three times as large, one can do:\n\njulia> C = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]\n3×3 Matrix{Float64}:\n 1.0  0.0  0.0\n 0.0  2.0  0.0\n 0.0  0.0  3.0\n\njulia> coupling = AnisoCoupling(C)\nAnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])\n\nOne can define effective 2D couplings by setting all coefficients in one of the dimensions to zero:\n\njulia> C = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]\n3x3 Matrix{Float64}\n 1.0 0.0 0.0\n 0.0 1.0 0.0\n 0.0 0.0 0.0\n\njulia> coupling_2d = AnisoCoupling(C)\nAnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0])\n\nSimilarly, one can define effective 1D couplings by setting two of the components to zero.\n\n\n\n\n\n","category":"type"},{"location":"couplingtensor/#SpiDy.IsoCoupling","page":"Spin-environment coupling tensor","title":"SpiDy.IsoCoupling","text":"struct IsoCoupling{TT<:Real} <: Coupling\n\nA subtype of Coupling used to represent isotropic coupling between spin and environment.\n\nFields\n\nC::Real: The coupling strength.\n\nExamples\n\njulia> coupling = IsoCoupling(2.5)\nIsoCoupling{Float64}(2.5)\n\n\n\n\n\n","category":"type"},{"location":"dynamics/#Dynamics","page":"Dynamics","title":"Dynamics","text":"","category":"section"},{"location":"dynamics/","page":"Dynamics","title":"Dynamics","text":"diffeqsolver","category":"page"},{"location":"dynamics/#SpiDy.diffeqsolver","page":"Dynamics","title":"SpiDy.diffeqsolver","text":"function diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], save_fields=false, projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)\n\nSolves the dynamics of a system of ineracting spins under the influence of both local (unique to each spin) and global (shared by all spins) stochastic noise from the environment.\n\nArguments\n\ns0: Array of length 3N specifying the initial conditions of the N spins. The order the initial consitions is first the Sx,Sy,Sz for the first spin, then for the second, and so on.\ntspan: The time span to solve the equations over, specified as a tuple (tstart, tend).\nJlist::Vector{LorentzianSD}: The list of spectral densities of the environment(s).\nbfield: An array of tuples of functions Array{Tuple{Function, Function, Function}} representing the time series of the noise for each environment.\nbcoupling::Vector{Vector{T}} T <: Real: A vector (of length the number of baths) of vectors (of length the number of spins), specifying if each bath couples to each spin.\nmatrix::Vector{TT} TT <: Coupling: A vector of Coupling structs specifying the spin-environment coupling matrix for each environment.\nJH=zero(I): (Optional) The spin-spin coupling matrix. Default is zero matrix (i.e. non-interacting spins).\nS0=1/2: (Optional) The spin quantum number. Default is 1/2.\nBext=[0, 0, 1]: (Optional) The external magnetic field vector. Default is [0, 0, 1] (normalised length pointing in the z direction).\nsaveat=[]: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.\nsave_fields=false: (Optional) If true, also return the auxiliary fields encoding the environment memory.\nprojection=false: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is false.\nalg=Tsit5(): (Optional) The differential equation solver algorithm. Default is Tsit5(). See the DifferentialEquations.jl docs for choices.\natol=1e-3: (Optional) The absolute tolerance for the solver. Default is 1e-3.\nrtol=1e-3: (Optional) The relative tolerance for the solver. Default is 1e-3.\n\nNote: The LorentzianSD type is provided by the SpectralDensities.jl package.\n\nNote: Additional keyword arguments will be passed on to the ODE solver (see the DifferentialEquations.jl docs)\n\nReturns\n\nAn ODESolution struct from DifferentialEquations.jl containing the solution of the equations of motion.\n\n\n\n\n\nfunction diffeqsolver(s0, tspan, J::LorentzianSD, Jshared::LorentzianSD, bfields, bfieldshared, matrix::Coupling; JH=zero(I), S0=1/2, Bext=[0, 0, 1], saveat=[], save_fields=false, projection=true, alg=Tsit5(), atol=1e-3, rtol=1e-3)\n\nSolves the dynamics of a system of ineracting spins under the influence of either local (unique to each spin) or global (shared by all spins) stochastic noise from the environment.\n\nArguments\n\ns0: Array of length 3N specifying the initial conditions of the N spins. The order the initial consitions is first the Sx,Sy,Sz for the first spin, then for the second, and so on.\ntspan: The time span to solve the equations over, specified as a tuple (tstart, tend).\nJ::LorentzianSD: The spectral density of the noise acting on the spins (either local or shared depending on the value of bfields).\nbfields: For local baths, an array of tuples of functions Array{Tuple{Function, Function, Function}} representing the time series of the local stochastic field for each spin. For a global bath, a tuple of functions Tuple{Function, Function, Function} representing the time series of the global stochastic field shared by all the spins.\nmatrix::Coupling: The spin-environment coupling matrix.\nJH=zero(I): (Optional) The spin-spin coupling matrix. Default is zero matrix (i.e. non-interacting spins).\nS0=1/2: (Optional) The spin quantum number. Default is 1/2.\nBext=[0, 0, 1]: (Optional) The external magnetic field vector. Default is [0, 0, 1] (normalised length pointing in the z direction).\nsaveat=[]: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.\nsave_fields=false: (Optional) If true, also return the auxiliary fields encoding the environment memory.\nprojection=false: (Optional) Specifies whether to project the spin vectors onto the unit sphere at each time step, hence forcing the numerical conservation of the spin length. Default is false.\nalg=Tsit5(): (Optional) The differential equation solver algorithm. Default is Tsit5(). See the DifferentialEquations.jl docs for choices.\natol=1e-3: (Optional) The absolute tolerance for the solver. Default is 1e-3.\nrtol=1e-3: (Optional) The relative tolerance for the solver. Default is 1e-3.\n\nNote: The LorentzianSD type is provided by the SpectralDensities.jl package.\n\nNote: Additional keyword arguments will be passed on to the ODE solver (see the DifferentialEquations.jl docs)\n\nReturns\n\nAn ODESolution struct from DifferentialEquations.jl containing the solution of the equations of motion.\n\n\n\n\n\nfunction diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfield, matrix::Coupling; JH=zero(I), Ω=1.0, counter_term=true, saveat=[], save_fields=false, alg=Tsit5(), atol=1e-3, rtol=1e-3)\n\nSolves the dynamics of a system of ineracting oscillators under the influence of either local (unique to each oscillator) or global (shared by all oscillators) stochastic noise from the environment.\n\nArguments\n\nx0: Array of length N specifying the initial position of the N oscillators.\np0: Array of length N specifying the initial momentum of the N oscillators.\ntspan: The time span to solve the equations over, specified as a tuple (tstart, tend).\nJ::LorentzianSD: The spectral density of the noise acting globally on the harmonic oscillators.\nbfield: A tuple of functions Tuple{Function, Function, Function} representing the time series of the global stochastic field shared by all the harmonic oscillators.\nmatrix::Coupling: The harmonic oscillators-environment coupling matrix.\nJH=zero(I): (Optional) The oscillator-oscillator coupling matrix. Default is zero matrix (i.e. non-interacting oscillators).\nΩ=1.0: (Optional) The natural angular frequency of the harmonic oscillators (currently the same for all). Default is 1.\ncounter_term=true: (Optional) Whether to include the counter-term or not. Default is true.\nsaveat=[]: (Optional) An array of time points where the solution should be saved. Default is empty, which saves the solution at the time steps chosen by the integration algorithm.\nsave_fields=false: (Optional) If true, also return the auxiliary fields encoding the environment memory.\nalg=Tsit5(): (Optional) The differential equation solver algorithm. Default is Tsit5(). See the DifferentialEquations.jl docs for choices.\natol=1e-3: (Optional) The absolute tolerance for the solver. Default is 1e-3.\nrtol=1e-3: (Optional) The relative tolerance for the solver. Default is 1e-3.\n\nNote: The LorentzianSD type is provided by the SpectralDensities.jl package.\n\nNote: Additional keyword arguments will be passed on to the ODE solver (see the DifferentialEquations.jl docs)\n\nReturns\n\nAn ODESolution struct from DifferentialEquations.jl containing the solution of the equations of motion.\n\n\n\n\n\n","category":"function"},{"location":"spectraldensity/#Spectral-density","page":"Spectral density","title":"Spectral density","text":"","category":"section"},{"location":"spectraldensity/","page":"Spectral density","title":"Spectral density","text":"Spectral densities are implemented using the SpectralDensities.jl package. Please refer to its documentation for details of how to use it. For SpiDy.jl, the key type is LorentzianSD, which is here used to characterise the spectral density of the environment.","category":"page"},{"location":"#SpiDy.jl-documentation","page":"Start with SpiDy.jl","title":"SpiDy.jl documentation","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"pages = [\"Start with SpiDy.jl\" => \"index.md\",\n         \"Units and parameter choice\" => \"units.md\",\n         \"Noise\" => \"noise.md\",\n         \"Spectral density\" => \"spectraldensity.md\",\n         \"Stochastic field\" => \"stochasticfield.md\",\n         \"Spin-environment coupling tensor\" => \"couplingtensor.md\",\n         \"Spin-spin interactions\" => \"couplingfunctions.md\",\n         \"Dynamics\" => \"dynamics.md]","category":"page"},{"location":"#Start-with-SpiDy.jl","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"SpiDy.jl is a Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper \"Quantum Brownian motion for magnets\" to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"The classical simulations in anisotropic coupling found in the pre-print \"Quantum-classical correspondence in spin-boson equilibrium states at arbitrary coupling\" have been generated using a very-early-version of this code.","category":"page"},{"location":"#Install-Julia","page":"Start with SpiDy.jl","title":"Install Julia","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are new to Julia, here is how to install it.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are a Windows/Mac user, download Julia here and run the installer. On Mac, drag-and-drop the app to the Applications.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are a Linux user, just open a terminal and use your package manager, e.g. on Debian-based distros run \"sudo apt-get install julia\", on RedHat-based distros run \"sudo dnf install julia\".","category":"page"},{"location":"#Install-SpiDy.jl","page":"Start with SpiDy.jl","title":"Install SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"Following the Julia General Registry guidelines, the package can be installed as follows. (NB: the entire installation of SpiDy and its dependencies takes about 5 minutes on a bare-bones Julia environment.)","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"Start Julia and enter in Pkg REPL mode by pressing ] then run the following,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"add https://github.com/quantum-exeter/SpiDy.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"or alternatively run the following lines in your code,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"using Pkg;\nPkg.add(url=\"https://github.com/quantum-exeter/SpiDy.jl\")","category":"page"},{"location":"#Run-SpiDy.jl","page":"Start with SpiDy.jl","title":"Run SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"To run the code,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"save run_dynamics.jl and run_steadystate.jl in your preferred location (right click -> save as... should work to save the file)\nopen the terminal or command line\nrun the following command,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"julia \"path-to-your-file\"/run_dynamics.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"where \"path-to-your-file\" is the one where you saved your file. Replace run_dynamics.jl with run_steadystate.jl to run the one of your choice.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"NB: the code can exploit parallel computation. To do this, run your files as","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"julia -t 6 \"path-to-your-file\"/run_dynamics.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"where you want to replace \"6\" with the number of threads that you wish to use. As a general idea, you do not want to use more than 80% of the number of threads you have available in your machine, e.g. if you have a 4-core CPU, you are likely to have 8 threads and you may want to run the parallelization as indicated above.","category":"page"},{"location":"#Units-and-choosing-parameters","page":"Start with SpiDy.jl","title":"Units and choosing parameters","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"See the section \"Units and choosing parameters\" for a discussion on the choice of units implemented in SpiDy.jl and how to choose these parameters given a concrete real problem.","category":"page"},{"location":"#Index","page":"Start with SpiDy.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"","category":"page"}]
}
