"""
    Nchain(N, J0; boundary=nothing)

Create the spin-spin coupling matrix of a 1D chain of `N` sites with nearest-neighbour interactions of strength set by `J0`.

# Arguments
- `N`: The number of sites in the chain.
- `J0`: The nearest-neighbour coupling strength.
- `boundary=nothing`: (Optional) Specifies the boundary condition of the chain. Default is `nothing` which corresponds to a chain open at the edges. Use `:periodic` for periodic boundary conditions.

# Returns
An `N×N` array representing the coupling matrix.

# Examples
```julia
julia> J = Nchain(4, 2.0)
4×4 Matrix{Float64}:
 0.0  2.0  0.0  0.0
 2.0  0.0  2.0  0.0
 0.0  2.0  0.0  2.0
 0.0  0.0  2.0  0.0

julia> J_periodic = Nchain(4, 2.0; boundary=:periodic)
4×4 Matrix{Float64}:
 0.0  2.0  0.0  2.0
 2.0  0.0  2.0  0.0
 0.0  2.0  0.0  2.0
 2.0  0.0  2.0  0.0
```
"""
function Nchain(N, J0; boundary=nothing)
    J = Array(SymTridiagonal(repeat([0.], N), repeat([J0], N-1)))
    if boundary == :periodic
        J[1,N] = J[N,1] = J0
    end
    return J
end

function map_1d_to_2d_idx(n, N)
    j, i = divrem(n-1, N)
    return i+1, j+1
end

"""
    NNlattice(N, Jh, Jv; boundary=nothing)

Create the spin-spin coupling matrix of a `NxN` 2D lattice with nearest-neighbour interactions and
specified horizontal and vertical coupling strengths, `Jh` and `Jv` respectively.

# Arguments
- `N`: The size of the lattice (N x N).
- `Jh`: The horizontal nearest-neighbour coupling strength.
- `Jv`: The vertical nearest-neighbour coupling strength.
- `boundary=nothing`: (Optional) Specifies the boundary condition of the lattice. Default is `nothing` which corresponds to open edges. Use `:periodic` for periodic boundary condition.

# Returns
An `N^2×N^2` array representing the coupling matrix of the lattice.

# Examples
```julia
julia> J = NNlattice(3, 2.0, 1.0)
9×9 LinearAlgebra.Symmetric{Float64, Matrix{Float64}}:
 0.0  2.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
 2.0  0.0  2.0  0.0  1.0  0.0  0.0  0.0  0.0
 0.0  2.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  2.0  0.0  1.0  0.0  0.0
 0.0  1.0  0.0  2.0  0.0  2.0  0.0  1.0  0.0
 0.0  0.0  1.0  0.0  2.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  1.0  0.0  0.0  0.0  2.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0  2.0  0.0  2.0
 0.0  0.0  0.0  0.0  0.0  1.0  0.0  2.0  0.0

julia> J_periodic = NNlattice(2, 2.0, 1.0; boundary=:periodic)
9×9 LinearAlgebra.Symmetric{Float64, Matrix{Float64}}:
 0.0  2.0  2.0  1.0  0.0  0.0  1.0  0.0  0.0
 2.0  0.0  2.0  0.0  1.0  0.0  0.0  1.0  0.0
 2.0  2.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0
 1.0  0.0  0.0  0.0  2.0  2.0  1.0  0.0  0.0
 0.0  1.0  0.0  2.0  0.0  2.0  0.0  1.0  0.0
 0.0  0.0  1.0  2.0  2.0  0.0  0.0  0.0  1.0
 1.0  0.0  0.0  1.0  0.0  0.0  0.0  2.0  2.0
 0.0  1.0  0.0  0.0  1.0  0.0  2.0  0.0  2.0
 0.0  0.0  1.0  0.0  0.0  1.0  2.0  2.0  0.0
```
"""
function NNlattice(N, Jh, Jv; boundary=nothing)
    J = zeros(N^2, N^2)
    for n in 1:N^2
        for m in 1:N^2
            jn, kn = map_1d_to_2d_idx(n, N)
            jm, km = map_1d_to_2d_idx(m, N)
            if (jn == jm && kn == km+1) || (jn == jm && kn == km-1) ||
               (boundary == :periodic && jn == jm && kn == 1 && km == N) ||
               (boundary == :periodic && jn == jm && kn == N && km == 1)
                J[n, m] = Jv
            elseif (jn == jm+1 && kn == km) || (jn == jm-1 && kn == km) ||
               (boundary == :periodic && jn == 1 && jm  == N && kn == km) ||
               (boundary == :periodic && jn == N && jm  == 1 && kn == km)
                J[n, m] = Jh
            end
        end
    end
    return Symmetric(J)
end
