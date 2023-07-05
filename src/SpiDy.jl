module SpiDy

using DifferentialEquations
using PreallocationTools
using Distributions
using FFTW
using Interpolations
using LinearAlgebra
using NPZ
using DataFrames
using CSV
using ProgressMeter
using Random
using Statistics

include("Noise.jl");
include("SpectralDensity.jl");
include("StochasticField.jl");
include("CouplingTensor.jl");
include("CouplingFunctions.jl");
include("Dynamics.jl");

# export external files to build documentation
export Noise, SpectralDensity, StochasticField, CouplingTensor, CouplingFunctions, Dynamics
# export structures to build documentation
export ClassicalNoise, QuantumNoise, NoZeroQuantumNoise, GenericSD, LorentzianSD, PolySD, Coupling, AnisoCoupling, IsoCoupling
# export modules and functions to build documentation
export spectrum, sd, sdoverω, reorgenergy, kernel, imagkernel, psd, bfield, Nchain, NNlattice, diffeqsolver

end
