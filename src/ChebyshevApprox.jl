function chebrow(n::Int)
    coeffs = zeros(Int, n)
    if n==1
        coeffs[1] += 1
        return coeffs
    elseif n==2
        coeffs[2] += 1
        return coeffs
    end
    coeffs[2:end] += 2*chebrow(n-1)
    coeffs[1:end-2] -= chebrow(n-2)
    return coeffs
end

function chebmatrix(n::Int)
    matrix = zeros(Int, (n, n))
    for i in 1:n
        matrix[i, :] = [chebrow(i); zeros(Int, n-i)]
    end
    return matrix
end

function polycoeffs(n, chebcoeffs)
    return chebmatrix(n)' * chebcoeffs
end