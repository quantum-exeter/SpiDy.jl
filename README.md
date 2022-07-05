# SpiDy.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://quantum-exeter.github.io/SpiDy.jl/dev)

Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper <a href=https://doi.org/10.1088/1367-2630/ac4ef2>"Quantum Brownian motion for magnets"</a> to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.

## Install the package
Following the Julia *General Registry*, the package can be installed as follows,
```Julia
>>> using Pkg;
>>> Pkg.add(url="https://github.com/quantum-exeter/SpiDy.jl")
```
or using the Pkg REPL mode as follows,
```Julia
>>> ]
>>> add https://github.com/quantum-exeter/SpiDy.jl
```

## Repo structure
* **.github/workflows**: contains the yml file to build the documentation and commit on the gh-pages branch
* **docs**: contains the logos, make.jl and index.md for the generation of documentation
* **notebooks**: contains jupyter notebooks to test the code
* **src**: contains the code

## Online documentation
Check the online documentation at <a href="https://quantum-exeter.github.io/SpiDy.jl/dev/">this link</a>.