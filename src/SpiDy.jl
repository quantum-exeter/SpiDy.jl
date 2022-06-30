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
export psd, b_field, spectrum, sd, sdoverω, reorgenergy, damping_kernel_frequency

end