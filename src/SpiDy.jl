module SpiDy

using LinearAlgebra
using DifferentialEquations
using FFTW
using Random

include("SpectralDensity.jl");
include("Noise.jl");
include("PSD.jl");

export SpectralDensity, LorentzianSD, Noise, PSD

end