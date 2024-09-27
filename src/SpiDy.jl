module SpiDy

using DifferentialEquations
using SpectralDensities
using PreallocationTools
using Distributions
using FFTW
using Interpolations
using LinearAlgebra
using Random
using Statistics

include("Noise.jl");
include("StochasticField.jl");
include("CouplingTensor.jl");
include("CouplingFunctions.jl");
include("Dynamics.jl");

# export external files to build documentation
export Noise, SpectralDensity, StochasticField, CouplingTensor, CouplingFunctions, Dynamics
# export structures to build documentation
export ClassicalNoise, QuantumNoise, NoZeroQuantumNoise, Coupling, AnisoCoupling, IsoCoupling
# export modules and functions to build documentation
export spectrum, psd, bfield, Nchain, NNlattice, diffeqsolver
# export useful SpectralDensities.jl types and functions
export SpectralDensities, AbstractSD, LorentzianSD

end
