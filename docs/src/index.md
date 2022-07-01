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
bfield(N, Δt, J::GenericSD, noise::Noise, distro=nothing)
```

<!-- ## Spin state
```@docs
```

## Coupling function
```@docs
```

## Dynamics
```@docs
diffeqsolver(N, Δt, J::GenericSD, noise::Noise, distro=nothing)
``` -->

## Index
```@index
```