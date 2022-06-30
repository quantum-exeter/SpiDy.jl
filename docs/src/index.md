# SpiDy.jl documentation
```@contents
```

## Noise

## Spectral density
```@docs
sd(J::GenericSD)
sdoverω(J::GenericSD)
sdoverω(J::LorentzianSD)
reorgenergy(J::GenericSD)
reorgenergy(J::LorentzianSD)
damping_kernel_frequency(J::LorentzianSD)
psd(J::GenericSD, noise::Noise)
psd(J::LorentzianSD, noise::ClassicalNoise)
```

## Stochastic field
```@docs
b_field(N, Δt, J::GenericSD, noise::Noise, distro=nothing)
```

## Index
```@index
```
