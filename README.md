![spidy_banner](https://github.com/quantum-exeter/SpiDy.jl/blob/main/docs/src/assets/banner.png)

# SpiDy.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantum-exeter.github.io/SpiDy.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://quantum-exeter.github.io/SpiDy.jl/dev/)
[![Build Status](https://github.com/quantum-exeter/SpiDy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/quantum-exeter/SpiDy.jl/actions/workflows/CI.yml?query=branch%3Amain)


Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper <a href=https://doi.org/10.1088/1367-2630/ac4ef2>"Quantum Brownian motion for magnets"</a> to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.

The classical simulations in anisotropic coupling found in the pre-print <a href=https://arxiv.org/abs/2204.10874>"Quantum-classical correspondence in spin-boson equilibrium states at arbitrary coupling"</a> have been generated using a very-early-version of this code.

## Install Julia
If you are new to Julia, here is how to install it.

If you are a Windows/Mac user, <a href=https://julialang.org/downloads/>download Julia here</a> and run the installer. On Mac, drag-and-drop the app to the Applications.

If you are a Linux user, just open a terminal and use your package manager, e.g. on Debian-based distros run "sudo apt-get install julia", on RedHat-based distros run "sudo dnf install julia".

## Install SpiDy
Following the Julia General Registry guidelines, the package can be installed as follows. *(NB: the entire installation of SpiDy and its dependencies takes about 5 minutes on a bare-bones Julia environment.)*

Start Julia and enter in Pkg REPL mode by pressing **]** then run the following,
```Julia
add https://github.com/quantum-exeter/SpiDy.jl
```
or alternatively run the following lines in your code,
```Julia
using Pkg;
Pkg.add(url="https://github.com/quantum-exeter/SpiDy.jl")
```

## Run SpiDy
To run the code,
* save <a href=https://raw.githubusercontent.com/quantum-exeter/SpiDy.jl/main/runs/run_dynamics.jl>run_dynamics.jl</a> and <a href=https://raw.githubusercontent.com/quantum-exeter/SpiDy.jl/main/runs/run_steadystate.jl>run_steadystate.jl</a> in your preferred location (open the link -> right click on the page -> save as... should work to save the file)
* open the terminal or command line
* run the following command,
```Julia
julia "path-to-your-file"/run_dynamics.jl
```
where "path-to-your-file" is the one where you saved your file. Replace *run_dynamics.jl* with *run_steadystate.jl* to run the one of your choice.

**This last command will run the code and save plots/datafile of the chosen run. CONGRATS, you have just run SpiDy for the first time!**

NB: the code can exploit parallel computation. To do this, run your files as
```Julia
julia -t 6 "path-to-your-file"/run_dynamics.jl
```
where you want to replace "6" with the number of threads that you wish to use. As a general idea, you do not want to use more than 80% of the number of threads you have available in your machine, e.g. if you have a 4-core CPU, you are likely to have 8 threads and you may want to run the parallelization as indicated above.


## Repo structure
* **.github/workflows**: contains the yml file to build the documentation and commit on the gh-pages branch
* **docs**: contains the logos, make.jl and index.md for the generation of documentation
* **runs**: contains run_*.jl files which can be used as a template to run the code
* **src**: contains the source code
* **starthere**: contains an ipynb notebook written in Julia which walks you through bits and pieces of the code with explanatory plots *(the notebook is evolving over time but always ready to use)*

## Online documentation
Check the online documentation at <a href="https://quantum-exeter.github.io/SpiDy.jl/dev/">this link</a>.
