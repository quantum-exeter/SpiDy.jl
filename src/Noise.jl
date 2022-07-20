"""
```Julia
Noise
```

Definition of the abstract type `Noise`.
"""
abstract type Noise end

"""
```Julia
ClassicalNoise{TT<:Real}
```

Returns a `ClassicalNoise` structure of type `Noise` built by passing a `Real` value.
This allows for future overloading of the `spectrum` method.
"""
struct ClassicalNoise{TT<:Real} <: Noise
    T::TT
end

"""
```Julia
ClassicalNoise{TT<:Real}
```

Returns a `QuantumNoise` structure of type `Noise` built by passing a `Real` value.
This allows for future overloading of the `spectrum` method.
"""
struct QuantumNoise{TT<:Real} <: Noise
    T::TT
end

"""
```Julia
ClassicalNoise{TT<:Real}
```

Returns a `NoZeroQuantumNoise` structure of type `Noise` built by passing a `Real` value.
This allows for future overloading of the `spectrum` method.
"""
struct NoZeroQuantumNoise{TT<:Real} <: Noise
    T::TT
end

"""
```Julia
spectrum(n::ClassicalNoise)
```

Returns the classical noise at temperature `n.T` as a function of `ω`. The conditional
statement makes sure the function does not divide by zero in case of `ω==0`.
"""
spectrum(n::ClassicalNoise) = ω -> iszero(ω) ? zero(ω) : 2*n.T/ω

"""
```Julia
spectrum(n::QuantumNoise)
```

Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
statement makes sure the function does not divide by zero in case of `ω==0`.
"""
spectrum(n::QuantumNoise) = ω -> iszero(n.T) ? sign(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2))

"""
```Julia
spectrum(n::NoZeroQuantumNoise)
```

Returns the quantum noise at temperature `n.T` as a function of `ω`. The conditional
statement makes sure the function does not divide by zero in case of `ω==0`. It differs
from `spectrum(n::QuantumNoise)` in the fact that the zero-point noise is removed.
"""
spectrum(n::NoZeroQuantumNoise) = ω -> iszero(n.T) ? zero(ω) : (iszero(ω) ? zero(ω) : coth(ω/n.T/2) - sign(ω))
