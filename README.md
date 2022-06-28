# SpiDy.jl
Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper <a href=https://doi.org/10.1088/1367-2630/ac4ef2>"Quantum Brownian motion for magnets"</a> to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.

## Repo structure
* **.github/workflows**: contains the yml file to build the documentation and commit on the gh-pages branch
* **docs**: contains the logos, make.jl and index.md for the generation of documentation
* **src**: contains the code

## Online documentation
Check the online documentation at <a href="https://quantum-exeter.github.io/SpiDy.jl/dev/">this link</a>.