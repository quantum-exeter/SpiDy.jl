push!(LOAD_PATH,"../src/")

using Documenter, SpiDy

makedocs(sitename="SpiDy.jl",
    pages = ["Index" => "index.md",
             "Noise" => "noise.md",
             "Spectral density" => "spectraldensity.md",
             "Stochastic field" => "stochasticfield.md",
             "Coupling tensor" => "couplingtensor.md",
             "Dynamics" => "dynamics.md"])

deploydocs(
    repo = "github.com/quantum-exeter/SpiDy.jl.git",
)