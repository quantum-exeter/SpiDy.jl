push!(LOAD_PATH,"../src/")

using Documenter, SpiDy

makedocs(sitename="SpiDy.jl")

deploydocs(
    repo = "github.com/quantum-exeter/SpiDy.jl.git",
)
