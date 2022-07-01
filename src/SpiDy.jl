module SpiDy

using LinearAlgebra
using DifferentialEquations
using FFTW
using Random
using Distributions

include("Noise.jl");
include("SpectralDensity.jl");
include("StochasticField.jl");

export SpectralDensity, Noise, StochasticField
export ClassicalNoise, QuantumNoise, NoZeroQuantumNoise, GenericSD, LorentzianSD
export psd, bfield, spectrum, sd, sdoverω, reorgenergy, kernel

end