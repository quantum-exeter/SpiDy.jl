## Spectral density
```@docs
GenericSD
LorentzianSD
PolySD
sd(J::GenericSD)
sdoverω(J::GenericSD)
sdoverω(J::LorentzianSD)
sd(J::PolySD)
reorgenergy(J::GenericSD)
reorgenergy(J::LorentzianSD)
kernel(J::LorentzianSD)
imagkernel(J::GenericSD)
psd(J::GenericSD, noise::Noise)
psd(J::LorentzianSD, noise::ClassicalNoise)
```