push!(LOAD_PATH,"../src/")

using Documenter, SpiDy

makedocs(sitename="SpiDy Documentation")

deploydocs(
    repo = "github.com/quantum-exeter/SpiDy.jl.git",
)
