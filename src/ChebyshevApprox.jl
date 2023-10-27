function cheb_row(n::Int)
    coeffs = zeros(Int, n)
    if n==1
        coeffs[1] += 1
        return coeffs
    elseif n==2
        coeffs[2] += 1
        return coeffs
    end
    coeffs[2:end] += 2*cheb_row(n-1)
    coeffs[1:end-2] -= cheb_row(n-2)
    return coeffs
end

function cheb_matrix(n::Int)
    matrix = zeros(Int, (n, n))
    for i in 1:n
        matrix[i, :] = [cheb_row(i); zeros(Int, n-i)]
    end
    return matrix
end

function poly_coeffs(n, cheb_coeffs)
    return cheb_matrix(n)' * cheb_coeffs
end