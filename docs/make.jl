push!(LOAD_PATH,"../src/")

using Documenter, SpiDy

makedocs(;
    modules=[SpiDy],
    authors=["Stefano Scali <scali.stefano@gmail.com>", "Federico Cerisola <federico@cerisola.net>"],
    repo="https://github.com/quantum-exeter/SpiDy.jl/blob/{commit}{path}#{line}",
    sitename="SpiDy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://quantum-exeter.github.io/SpiDy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages = ["Start with SpiDy" => "startwithspidy.md",
             "Index" => "index.md",
             "Noise" => "noise.md",
             "Spectral density" => "spectraldensity.md",
             "Stochastic field" => "stochasticfield.md",
             "Coupling tensor" => "couplingtensor.md",
             "Coupling functions" => "couplingfunctions.md",
             "Dynamics" => "dynamics.md"],
)

deploydocs(
    repo = "github.com/quantum-exeter/SpiDy.jl.git",
    devbranch="main",
)