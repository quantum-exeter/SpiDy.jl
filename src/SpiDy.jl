module SpiDy

using LinearAlgebra
using DifferentialEquations
using FFTW
using Random
using Distributions

include("SpectralDensity.jl");
include("Noise.jl");
include("StochasticField.jl");

export SpectralDensity, LorentzianSD, Noise, PSD

end