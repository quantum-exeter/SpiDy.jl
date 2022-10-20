"""
```Julia
NNChain(N, J0=1; boundary=nothing)
```

Returns the coupling matrix for a 1D chain of `N` sites with nearest-neighbours interactions of strength set by `J0`.

Keyword arguments:
- `boundary` if set to `:periodic` then periodic boundary conditions are considered, otherwise the chain is open at the edges
```
"""
function NNChain(N, J0=1; boundary=nothing)
    if boundary == :periodic
        J = zeros(N, N)
        for i in 1:N
            for j in 1:N
                if (i == j+1) || (i == j-1) || (i == 1 && j == N) || (i == N && j == 1)
                    J[i,j] = J0
                end
            end
        end
        return Symmetric(J)
    else
        return SymTridiagonal(repeat([0], N), repeat([J0], N-1))
    end
end

function map_1d_to_2d_idx(n, L)
    j, i = divrem(n-1, L)
    return i+1, j+1
end

"""
```Julia
NNSquareLattice(L, Jh=1, Jv=1; boundary=nothing)
```

Returns the coupling matrix for a 2D square lattice of `L`x`L` sites with nearest-neighbours interactions of strength `Jh` in the horizontal direction and `Jv` in the vertical one.

Keyword arguments:
- `boundary` if set to `:periodic` then periodic boundary conditions are considered, otherwise the lattice is open at the edges
```
"""
function NNSquareLattice(L, Jh=1, Jv=1; boundary=nothing)
    J = zeros(L^2, L^2)
    for n in 1:L^2
        for m in 1:L^2
            jn, kn = map_1d_to_2d_idx(n, L)
            jm, km = map_1d_to_2d_idx(m, L)
            if (jn == jm && kn == km+1) || (jn == jm && kn == km-1) ||
               (boundary == :periodic && jn == jm && kn == 1 && km == L) ||
               (boundary == :periodic && jn == jm && kn == L && km == 1)
                J[n, m] = Jv
            elseif (jn == jm+1 && kn == km) || (jn == jm-1 && kn == km) ||
               (boundary == :periodic && jn == 1 && jm  == L && kn == km) ||
               (boundary == :periodic && jn == L && jm  == 1 && kn == km)
                J[n, m] = Jh
            end
        end
    end
    return Symmetric(J)
end
