var documenterSearchIndex = {"docs":
[{"location":"couplingfunctions/#Spin-spin-interactions","page":"Coupling functions","title":"Spin-spin interactions","text":"","category":"section"},{"location":"couplingfunctions/","page":"Coupling functions","title":"Coupling functions","text":"Nchain(N, J0=1; boundary=nothing)\nNNlattice(L, Jh=1, Jv=1; boundary=nothing)","category":"page"},{"location":"noise/#Noise","page":"Noise","title":"Noise","text":"","category":"section"},{"location":"noise/","page":"Noise","title":"Noise","text":"Noise\nClassicalNoise\nQuantumNoise\nNoZeroQuantumNoise\nspectrum(n::ClassicalNoise)\nspectrum(n::QuantumNoise)\nspectrum(n::NoZeroQuantumNoise)","category":"page"},{"location":"noise/#SpiDy.Noise","page":"Noise","title":"SpiDy.Noise","text":"abstract type Noise\n\nAn abstract type used to represent different kinds of environment noise.\n\nAny user-defined subtype of Noise should implement a method noise(::NoiseSubType) which should return a function that represent the spectrum of the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.ClassicalNoise","page":"Noise","title":"SpiDy.ClassicalNoise","text":"struct ClassicalNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent classical noise, with spectrum\n\nmathcalS(omega) = frack_mathrmBThbaromega\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.QuantumNoise","page":"Noise","title":"SpiDy.QuantumNoise","text":"struct QuantumNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent quantum noise, with spectrum\n\nmathcalS(omega) = cothleft(frachbaromega2k_mathrmBTright)\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"noise/#SpiDy.NoZeroQuantumNoise","page":"Noise","title":"SpiDy.NoZeroQuantumNoise","text":"struct NoZeroQuantumNoise{TT<:Real} <: Noise\n\nA subtype of Noise used to represent quantum noise with zero-point-fluctuations removed, that is with spectrum\n\nmathcalS(omega) = cothleft(frachbaromega2k_mathrmBTright) - 1\n\nFields\n\nT::TT: The temperature of the environment associated with the noise.\n\n\n\n\n\n","category":"type"},{"location":"stochasticfield/#Stochastic-field","page":"Stochastic field","title":"Stochastic field","text":"","category":"section"},{"location":"stochasticfield/","page":"Stochastic field","title":"Stochastic field","text":"bfield(N, Δt, J::GenericSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true)","category":"page"},{"location":"couplingtensor/#Environment-coupling-tensor","page":"Coupling tensor","title":"Environment coupling tensor","text":"","category":"section"},{"location":"couplingtensor/","page":"Coupling tensor","title":"Coupling tensor","text":"Coupling\nAnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}}\nIsoCoupling{TT<:Real}","category":"page"},{"location":"couplingtensor/#SpiDy.Coupling","page":"Coupling tensor","title":"SpiDy.Coupling","text":"abstract type Coupling\n\nAbstract type used to represent different forms of the environment coupling tensor.\n\n\n\n\n\n","category":"type"},{"location":"couplingtensor/#SpiDy.AnisoCoupling","page":"Coupling tensor","title":"SpiDy.AnisoCoupling","text":"struct AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}} <: Coupling\n\nA subtype of Coupling used to represent generic anisotropic coupling between the spin and environment.\n\nFields\n\nC::TT: The coupling matrix, which must be a 3x3 matrix (of type Matrix{<:Real}).\n\nExamples\n\nFor a system with environment coupling where the y coupling is twice as large as the x, and the z three times as large, one can do:\n\njulia> C = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]\n3×3 Matrix{Float64}:\n 1.0  0.0  0.0\n 0.0  2.0  0.0\n 0.0  0.0  3.0\n\njulia> coupling = AnisoCoupling(C)\nAnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])\n\nOne can define effective 2D couplings by setting all coefficients in one of the dimensions to zero:\n\njulia julia> C = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0] 3x3 Matrix{Float64}  1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 0.0\n\njulia> coupling_2d = AnisoCoupling(C) AnisoCoupling{Matrix{Float64}}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]) ``` Similarly, one can define effective 1D couplings by setting two of the components to zero.\n\n\n\n\n\n","category":"type"},{"location":"couplingtensor/#SpiDy.IsoCoupling","page":"Coupling tensor","title":"SpiDy.IsoCoupling","text":"struct IsoCoupling{TT<:Real} <: Coupling\n\nA subtype of Coupling used to represent isotropic coupling between spin and environment.\n\nFields\n\nC::Real: The coupling strength.\n\nExamples\n\njulia> coupling = IsoCoupling(2.5)\nIsoCoupling{Float64}(2.5)\n\n\n\n\n\n","category":"type"},{"location":"dynamics/#Dynamics","page":"Dynamics","title":"Dynamics","text":"","category":"section"},{"location":"dynamics/","page":"Dynamics","title":"Dynamics","text":"diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1], saveat=[])\ndiffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])","category":"page"},{"location":"dynamics/#SpiDy.diffeqsolver-Tuple{Any, Any, Any, LorentzianSD, Any, Coupling}","page":"Dynamics","title":"SpiDy.diffeqsolver","text":"diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; JH=zero(I), Ω=1.0, saveat=[], alg=Tsit5(), atol=1e-3, rtol=1e-3)\n\nReturns [sol.t, x, p, xinterp, pinterp], that is,\n\nthe vector sol.t of time steps at which the solutions are evaluated,\nthe vectors of the solutions positions x and momenta p evaluated at times sol.t,\nthe functions xinterp(t), pinterp(t) interpolations of the solutions found in the given time span.\n\nThe differential equation solver is built to account for Lorentzian spectral density.\n\nKeyword arguments:\n\nΩ harmonic oscillator bare frequency set as default to Ω=1.0.\nsaveat is an option of the function solve() which allows to only save the solution at the points needed to evaluate the steady-state, i.e. at late times. Used to optimize memory management and speed of the solver. Default value is an empty list, saveat=[], resulting in the solution being saved at optimal time steps within the time span.\nJH is the Heisenberg coupling matrix. Note that this have to be a symmetric matrix with zero diagonal. The preset value is the additive identity of the UniformScaling type, JH=zero(I).\nalg chooses the solving algorithm.\natol and rtol define the absolute and relative tolerance of the solver.\n\nExamples\n\njulia> diffeqsolver(x0, p0, tspan, J, bfields, matrix; saveat=((N*4÷5):1:N)*Δt)\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#Spectral-density","page":"Spectral density","title":"Spectral density","text":"","category":"section"},{"location":"spectraldensity/","page":"Spectral density","title":"Spectral density","text":"GenericSD\nLorentzianSD\nPolySD\nsd(J::GenericSD)\nsdoverω(J::GenericSD)\nsdoverω(J::LorentzianSD)\nsd(J::PolySD)\nreorgenergy(J::GenericSD)\nreorgenergy(J::LorentzianSD)\nkernel(J::LorentzianSD)\nimagkernel(J::GenericSD)\npsd(J::GenericSD, noise::Noise)\npsd(J::LorentzianSD, noise::ClassicalNoise)","category":"page"},{"location":"spectraldensity/#SpiDy.GenericSD","page":"Spectral density","title":"SpiDy.GenericSD","text":"GenericSD\n\nDefinition of the abstract type GenericSD.\n\n\n\n\n\n","category":"type"},{"location":"spectraldensity/#SpiDy.LorentzianSD","page":"Spectral density","title":"SpiDy.LorentzianSD","text":"LorentzianSD{T<:Real}\n\nReturns a LorentzianSD structure of type GenericSD built by passing 3 Real values. The values are ordered as LorentzianSD(α, ω0, Γ).\n\nExamples\n\njulia> LorentzianSD(1., 3., 8.)\n\n\n\n\n\n","category":"type"},{"location":"spectraldensity/#SpiDy.PolySD","page":"Spectral density","title":"SpiDy.PolySD","text":"PolySD{T<:Real}\n\nReturns a PolySD structure of type GenericSD built by passing 3 Real values. The values are ordered as PolySD(s, α, ωcut).\n\nExamples\n\njulia> PolySD(1., 3., 100.)\n\n\n\n\n\n","category":"type"},{"location":"spectraldensity/#SpiDy.sd-Tuple{GenericSD}","page":"Spectral density","title":"SpiDy.sd","text":"sd(J::GenericSD)\n\nDefines the spectral density for generic shapes GenericSD. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.sdoverω-Tuple{GenericSD}","page":"Spectral density","title":"SpiDy.sdoverω","text":"sdoverω(J::GenericSD)\n\nReturns the spectral density divided by ω for generic shapes GenericSD. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.sdoverω-Tuple{LorentzianSD}","page":"Spectral density","title":"SpiDy.sdoverω","text":"sdoverω(J::LorentzianSD)\n\nReturns the spectral density divided by ω for LorentzianSD shapes which naturally defines sd(J::LorentzianSD). The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.sd-Tuple{PolySD}","page":"Spectral density","title":"SpiDy.sd","text":"sdoverω(J::PolySD)\n\nReturns the spectral density for PolySD shapes with exponential cut-off which naturally defines sdoverω(J::PolySD). The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.reorgenergy-Tuple{GenericSD}","page":"Spectral density","title":"SpiDy.reorgenergy","text":"reorgenergy(J::GenericSD)\n\nReturns the reorganization energy numerically integrated as int_0^infty textsdoverω(omega)domega.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.reorgenergy-Tuple{LorentzianSD}","page":"Spectral density","title":"SpiDy.reorgenergy","text":"reorgenergy(J::LorentzianSD)\n\nReturns the analytical reorganization energy for LorentzianSD shapes.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.kernel-Tuple{LorentzianSD}","page":"Spectral density","title":"SpiDy.kernel","text":"kernel(J::LorentzianSD)\n\nReturns the specific damping kernel for a Lorentzian spectral density defined by the parameters in J. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.imagkernel-Tuple{GenericSD}","page":"Spectral density","title":"SpiDy.imagkernel","text":"imagkernel(J::GenericSD)\n\nReturns the imaginary part of the damping kernel for a Generic spectral density real and anti-symmetric in ω. Of this kind, we find Lorentzian spectral densities, and Polynomial spectral densities with sinmathbbN. The spectral density is defined by the parameters in J. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.psd-Tuple{GenericSD, Noise}","page":"Spectral density","title":"SpiDy.psd","text":"psd(J::GenericSD, noise::Noise)\n\nReturns the power spectral density depending on parameters J and noise. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"spectraldensity/#SpiDy.psd-Tuple{LorentzianSD, ClassicalNoise}","page":"Spectral density","title":"SpiDy.psd","text":"psd(J::LorentzianSD, noise::ClassicalNoise)\n\nReturns the analytical expression for power spectrum depending on Lorentzian spectral density and Classical noise. The returned function depends on ω.\n\n\n\n\n\n","category":"method"},{"location":"#SpiDy.jl-documentation","page":"Start with SpiDy.jl","title":"SpiDy.jl documentation","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"pages = [\"Start with SpiDy.jl\" => \"index.md\",\n         \"Noise\" => \"noise.md\",\n         \"Spectral density\" => \"spectraldensity.md\",\n         \"Stochastic field\" => \"stochasticfield.md\",\n         \"Coupling tensor\" => \"couplingtensor.md\",\n         \"Dynamics\" => \"dynamics.md\"])","category":"page"},{"location":"#Start-with-SpiDy.jl","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"SpiDy.jl is a Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper \"Quantum Brownian motion for magnets\" to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"The classical simulations in anisotropic coupling found in the pre-print \"Quantum-classical correspondence in spin-boson equilibrium states at arbitrary coupling\" have been generated using a very-early-version of this code.","category":"page"},{"location":"#Install-Julia","page":"Start with SpiDy.jl","title":"Install Julia","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are new to Julia, here is how to install it.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are a Windows/Mac user, download Julia here and run the installer. On Mac, drag-and-drop the app to the Applications.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"If you are a Linux user, just open a terminal and use your package manager, e.g. on Debian-based distros run \"sudo apt-get install julia\", on RedHat-based distros run \"sudo dnf install julia\".","category":"page"},{"location":"#Install-SpiDy.jl","page":"Start with SpiDy.jl","title":"Install SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"Following the Julia General Registry guidelines, the package can be installed as follows. (NB: the entire installation of SpiDy and its dependencies takes about 5 minutes on a bare-bones Julia environment.)","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"Start Julia and enter in Pkg REPL mode by pressing ] then run the following,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"add https://github.com/quantum-exeter/SpiDy.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"or alternatively run the following lines in your code,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"using Pkg;\nPkg.add(url=\"https://github.com/quantum-exeter/SpiDy.jl\")","category":"page"},{"location":"#Run-SpiDy.jl","page":"Start with SpiDy.jl","title":"Run SpiDy.jl","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"To run the code,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"save run_dynamics.jl and run_steadystate.jl in your preferred location (right click -> save as... should work to save the file)\nopen the terminal or command line\nrun the following command,","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"julia \"path-to-your-file\"/run_dynamics.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"where \"path-to-your-file\" is the one where you saved your file. Replace run_dynamics.jl with run_steadystate.jl to run the one of your choice.","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"NB: the code can exploit parallel computation. To do this, run your files as","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"julia -t 6 \"path-to-your-file\"/run_dynamics.jl","category":"page"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"where you want to replace \"6\" with the number of threads that you wish to use. As a general idea, you do not want to use more than 80% of the number of threads you have available in your machine, e.g. if you have a 4-core CPU, you are likely to have 8 threads and you may want to run the parallelization as indicated above.","category":"page"},{"location":"#Index","page":"Start with SpiDy.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Start with SpiDy.jl","title":"Start with SpiDy.jl","text":"","category":"page"}]
}
