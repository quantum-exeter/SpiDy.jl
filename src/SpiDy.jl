module SpiDy

using LinearAlgebra
using DifferentialEquations
using FFTW
using Random
using Distributions

include("SpectralDensity.jl");
include("Noise.jl");
include("StochasticField.jl");

export SpectralDensity, Noise, StochasticField
export LorentzianSD, PSD, b_field

end