## Spectral density
```@docs
GenericSD
LorentzianSD
sd(J::GenericSD)
sdoverω(J::GenericSD)
sdoverω(J::LorentzianSD)
reorgenergy(J::GenericSD)
reorgenergy(J::LorentzianSD)
kernel(J::LorentzianSD)
psd(J::GenericSD, noise::Noise)
psd(J::LorentzianSD, noise::ClassicalNoise)
```