## Units and choice of input parameters

In this document we describe the choice of units implemented in SpiDy.jl (which
follows the conventions utilised in 
**[NJP 24 033020 (2022)](https://www.doi.org/10.1088/1367-2630/ac4ef2)**) and
how to appropriately choose the parameters that the library takes as input.

### Equations of motion

SpiDy.jl is designed to implement the equations of spins in presence of a bath
of harmonic oscillators with Lorentzian spectral density, as derived in
**[NJP 24 033020 (2022)](https://www.doi.org/10.1088/1367-2630/ac4ef2)**, which
*for a single spin ``\mathbf{S}`` in presence of an external field
``\mathbf{B}_\mathrm{ext}`` read
```math
\frac{\mathrm{d}\mathbf{S}}{\mathrm{d}t} =
    \gamma\mathbf{S}\times\left(\mathbf{S}_\mathrm{ext} + \mathbf{b} + \mathbf{V}\right), \\
\frac{\mathrm{d}\mathbf{V}}{\mathrm{d}t} = \mathbf{W}, \\
\frac{\mathrm{d}\mathbf{W}}{\mathrm{d}t} = \gamma A \mathbf{S} - \omega_0^2\mathbf{V} - \Gamma\mathbf{W},
```
where ``A``, ``\omega_0``, and ``\Gamma`` parametrise the Lorentzian spectral density as
```math
J(\omega) = \frac{A\Gamma}{\pi} \frac{\omega}{(\omega_0^2 - \omega^2)^2 + \omega^2\Gamma^2},
```
and the stochastic ``\mathbf{b}`` is given by
```math
\mathbf{b}(t) = \int_{-\infty}^{+\infty}\mathrm{d}t' F(t-t') \xi(t'),
```
with ``\xi`` being white noise and
```math
F(\tau) = \frac{1}{2\pi}\int_{-\infty}^{+\infty}\mathrm{d}\omega
    e^{-i\omega\tau} \sqrt{P(\omega)}.
```
Here, ``P(\omega)`` is the power spectral density of the environment and it is
given in terms of the Lorentzian spectral density ``J(\omega)`` and the
environment thermal noise ``N(\omega)`` by ``P(\omega) = \hbar\pi J(\omega)
N(\omega)``. The noise ``N(\omega)`` can be classical, `quantum`, or `quantum` with
no zero point fluctuations. For example, for the `quantum` case we have
```math
N_\mathrm{qu}(\omega) = \coth\left(\frac{\hbar\omega}{2k_\mathrm{B}T}\right).
```

Note that in these equations above, all quantities have standard units.

### Units rescaling

SpiDy.jl implements these equations in a unit-free way, following the
conventions of **[NJP 24 033020 (2022)](https://www.doi.org/10.1088/1367-2630/ac4ef2)**.
In summary, suppose a spin of length ``\hbar S_0`` is in presence of a magnetic
field of magnitude ``B_0``. We define the Larmor frequency
```math
\omega_\mathrm{L} = |\gamma_e| B_0,
```
with ``\gamma_e`` the electron gyromagnetic ratio.

SpiDy.jl then takes as input parameters for the simulation the following **unit free**
parameters:
- ``S_0``: the spin length (in units of ``\hbar``).
- ``\bar{B}_\mathrm{ext}``: external magnetic field.
- ``\bar{t}_\mathrm{end}``: final time of the evolution.
- ``\mathrm{d}\bar{t}``: time differential.
- ``\bar{\omega}_0``: peak frequency of the Lorentzian spectral density.
- ``\bar{\Gamma}``: with of the Lorentzian spectral density.
- ``\bar{\alpha}``: amplitude of the Lorentzian density.
- ``\bar{T}``: the environment temperature.

These quantities are related to the unitful units in the previous section by the
following conversions:
```math
\begin{align}
\mathbf{B}_\mathrm{ext} &= B_0 \, \bar{B}_\mathrm{ext}, \\
t_\mathrm{end} &= \omega_\mathrm{L}^{-1} \, \bar{t}_\mathrm{end}, \\
\mathrm{d}t &= \omega_\mathrm{L}^{-1} \, \mathrm{d}\bar{t}, \\
\omega_0 &= \omega_\mathrm{L} \, \bar{\omega}_0, \\
\Gamma &= \omega_\mathrm{L} \, \bar{\Gamma}, \\
\omega_0 &= \omega_\mathrm{L} \, \bar{\omega}_0, \\
A &= \frac{B_0^2\omega_\mathrm{L}}{\hbar S_0} \, \bar{\alpha}, \\
T &= \frac{\hbar\omega_\mathrm{L}}{k_\mathrm{B}} \, \bar{T}.
\end{align}
```

With these definitions, the unit-free Gilbert damping is given by
```math
\eta = \frac{\bar{\alpha}\bar{\Gamma}}{\bar{\omega}_0^4}.
```

!!! warning "Units convetion in old versions of SpiDy.jl"
    Note that the unit conversion just explained, in line with
    **[NJP 24 033020 (2022)](https://www.doi.org/10.1088/1367-2630/ac4ef2)**,
    is correct for SpiDy.jl versions 1.2.0 and up.

    Older version of SpiDy.jl used a different convention, which unfortunately
    ment that a change of spin length ``S_0`` implied a change in temperature and
    time-scale, therefore requiring redefinig the paramerts of the Loretnzian,
    evolution time, etc.
    
    The change in units convetion means that versions of SpiDy.jl before and
    after 1.2.0 will produce different results. This breaking change was made to
    bring the code into consistency with the article
    **[NJP 24 033020 (2022)](https://www.doi.org/10.1088/1367-2630/ac4ef2)** and
    to make the mapping of parameters used in SpiDy.jl to real units much easier
    and more straightforward (see also next section).

    Finally, it is worth noting that the behaviour of SpiDy.jl in **older
    versions** can be exactly **recovered** by passing the option **``S_0 = 1``**
    to the integrator, since in that case both units **conventions agree**.

## Extracting Lorentzian units from experimental data

To be written.