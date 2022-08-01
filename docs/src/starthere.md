# Start with SpiDy
Spin-Dynamics Julia package. The code is a generalization of the results obtained in the paper ["Quantum Brownian motion for magnets"](https://doi.org/10.1088/1367-2630/ac4ef2) to account for arbitrary dimensional system-bath coupling. The system considered is a quantized three-dimensional spin + environment Hamiltonian. The code solves a set of differential equations for the spin vector where the damping accounts for memory, arbitrary noise and arbitrary statistics.

The classical simulations in anisotropic coupling found in the pre-print ["Quantum-classical correspondence in spin-boson equilibrium states at arbitrary coupling"](https://arxiv.org/abs/2204.10874) have been generated using a very-early-version of this code.

## Install Julia
If you are new to Julia, here is how to install it.

If you are a Windows/Mac user, [download Julia here](https://julialang.org/downloads/) and run the installer. On Mac, drag-and-drop the app to the Applications.

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
* save [run_dynamics.jl](https://raw.githubusercontent.com/quantum-exeter/SpiDy.jl/main/runs/run_dynamics.jl) and [run_steadystate.jl](https://raw.githubusercontent.com/quantum-exeter/SpiDy.jl/main/runs/run_steadystate.jl) in your preferred location (right click -> save as... should work to save the file)
* open the terminal or command line
* run the following command,
```Julia
julia "path-to-your-file"/run_dynamics.jl
```
where "path-to-your-file" is the one where you saved your file. Replace *run_dynamics.jl* with *run_steadystate.jl* to run the one of your choice.

NB: the code can exploit parallel computation. To do this, run your files as
```Julia
julia -t 6 "path-to-your-file"/run_dynamics.jl
```
where you want to replace "6" with the number of threads that you wish to use. As a general idea, you do not want to use more than 80% of the number of threads you have available in your machine, e.g. if you have a 4-core CPU, you are likely to have 8 threads and you may want to run the parallelization as indicated above.