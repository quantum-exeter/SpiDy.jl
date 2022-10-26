"""
```Julia
Nchain(N, J0; boundary=nothing)
```

Returns the coupling matrix for a 1D chain of `N` sites with nearest-neighbours interactions of strength set by `J0`.

Keyword arguments:
- `boundary` if set to `:periodic` then periodic boundary conditions are considered, otherwise the chain is open at the edges
```
"""
function Nchain(N, J0; boundary=nothing)
    J = SymTridiagonal(repeat([0.], N), repeat([J0], N-1))
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
```Julia
NNlattice(N, Jh, Jv; boundary=nothing)
```

Returns the coupling matrix for a 2D square lattice of `N`x`N` sites with nearest-neighbours interactions of strength `Jh` in the horizontal direction and `Jv` in the vertical one.

Keyword arguments:
- `boundary` if set to `:periodic` then periodic boundary conditions are considered, otherwise the lattice is open at the edges
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
