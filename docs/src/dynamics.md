## Dynamics
```@docs
diffeqsolver(s0, tspan, J::LorentzianSD, bfields, matrix::Coupling; S0=1/2, Bext=[0, 0, 1], saveat=[])
diffeqsolver(x0, p0, tspan, J::LorentzianSD, bfields, matrix::Coupling; Ω=1.0, saveat=[])
```