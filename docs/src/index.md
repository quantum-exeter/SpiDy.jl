# SpiDy.jl documentation
```@contents
```

## Noise
```@docs
spectrum(n::ClassicalNoise)
spectrum(n::QuantumNoise)
spectrum(n::NoZeroQuantumNoise)
```

## Spectral density
```@docs
sd(J::GenericSD)
sdoverω(J::GenericSD)
sdoverω(J::LorentzianSD)
reorgenergy(J::GenericSD)
reorgenergy(J::LorentzianSD)
kernel(J::LorentzianSD)
psd(J::GenericSD, noise::Noise)
psd(J::LorentzianSD, noise::ClassicalNoise)
```

## Stochastic field
```@docs
bfield(N, Δt, J::GenericSD, noise::Noise; distro=Normal(0., 1/sqrt(Δt)), interpolation=true)
```

## Coupling tensor
```@docs
AnisoCoupling{TT<:AbstractMatrix{T} where {T<:Real}}
IsoCoupling{TT<:Real}
```

## Dynamics
```@docs
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1])
```

## Index
```@index
```