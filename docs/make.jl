push!(LOAD_PATH,"../src/")

using Documenter, SpiDy

makedocs(sitename="SpiDy.jl",
    pages = ["Index" => "src/index.md",
             "Noise" => "src/noise.md",
             "Spectral density" => "src/spectraldensity.md",
             "Stochastic field" => "src/stochasticfield.md",
             "Coupling tensor" => "src/couplingtensor.md",
             "Dynamics" => "src/dynamics.md"])

deploydocs(
    repo = "github.com/quantum-exeter/SpiDy.jl.git",
)