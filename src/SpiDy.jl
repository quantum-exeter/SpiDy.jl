module SpiDy

using DifferentialEquations
using Distributions
using FFTW
using Interpolations
using LinearAlgebra
using NPZ
using DataFrames
using CSV
using Plots
using ProgressMeter
using Random
using Statistics

include("Noise.jl");
include("SpectralDensity.jl");
include("StochasticField.jl");
include("CouplingTensor.jl");
include("Dynamics.jl");

# export external files to build documentation
export Noise, SpectralDensity, StochasticField, CouplingTensor, Dynamics
# export structures to build documentation
export ClassicalNoise, QuantumNoise, NoZeroQuantumNoise, GenericSD, LorentzianSD, Coupling, AnisoCoupling, IsoCoupling
# export modules and functions to build documentation
export psd, bfield, spectrum, sd, sdoverω, reorgenergy, kernel, diffeqsolver

end