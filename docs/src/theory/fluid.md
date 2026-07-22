# Fluid convective resistance

The fluid convective resistance ``R_f`` accounts for the boundary layer between the moving
fluid and the inner pipe wall. Computing it requires the convective heat transfer coefficient,
which in turn follows from the flow regime through the Reynolds, Prandtl and Nusselt numbers.
This page covers the full chain, all implemented in
[`resistance_fluid.jl`](https://github.com/GHE-jl/BoreholeResistance.jl/blob/main/src/resistance_fluid.jl).

## Reynolds number

The Reynolds number classifies the flow regime. For a circular pipe the characteristic length
is the hydraulic diameter ``D = 2r``:

```math
Re = \frac{\rho_f \, \dot V \, D}{\mu_f} = \frac{2 r \, \rho_f \, \dot V}{\mu_f},
```

where ``\dot{V}`` is the **mean fluid speed** in m/s (not the volumetric flow rate ``V``), ``r`` the
pipe inner radius, ``\rho_f`` the density and ``\mu_f`` the dynamic viscosity. The same formula
applies to an annular pipe system by using its hydraulic radius ``r = r_b - r_o``. See
[`Reynolds`](@ref).

## Prandtl number

The Prandtl number is the ratio of momentum to thermal diffusivity:

```math
Pr = \frac{c_f \, \mu_f}{k_f}.
```

[`Prandtl`](@ref) provides two overloads: one taking the mass-specific heat ``c_f`` [J/kg·K]
directly, and one taking the volumetric specific heat ``C_f = c_f \rho_f`` [J/m³·K], in which
case ``Pr = C_f \mu_f / (k_f \rho_f)``. Both evaluate to the same number.

## Friction factor

The Darcy friction factor ``f`` enters the turbulent Nusselt correlation. In the **laminar**
regime it is exact:

```math
f = \frac{64}{Re}, \qquad Re < 2300.
```

In the **turbulent** regime the package offers two routes.

### Colebrook–White (implicit)

The standard implicit correlation, solved by fixed-point iteration:

```math
\frac{1}{\sqrt{f}} = -2 \log_{10}\!\left(\frac{\epsilon}{3.7\,D} + \frac{2.51}{Re\,\sqrt{f}}\right),
```

with ``\epsilon`` the pipe roughness and ``D = 2r``. Implemented in
[`friction_factor_Colebrook_White`](@ref); the iteration converges when successive iterates
agree to ``10^{-5}``.

### Tkachenko–Mileikovskyi (explicit)

An explicit approximation that avoids iteration and agrees with Colebrook–White to within ~1 %
under typical GHE flow conditions — convenient where many evaluations are needed. See
[`friction_factor_Tkachenko_Mileikovskyi`](@ref).

## Nusselt number — Gnielinski correlation

The Nusselt number ``Nu = hD / k_f`` gives the convective coefficient ``h``. The package uses
the Gnielinski framework with explicit handling of the three regimes (Lamarche, 2023):

```math
Nu =
\begin{cases}
4 & Re < 2300 \quad\text{(laminar)} \\[6pt]
(1-\gamma)\,4 + \gamma\,Nu_{4000} & 2300 \le Re < 4000 \quad\text{(transition)} \\[6pt]
\dfrac{(f/8)\,(Re-1000)\,Pr}{1 + 12.7\,\sqrt{f/8}\,\bigl(Pr^{2/3}-1\bigr)} & Re \ge 4000 \quad\text{(turbulent)}
\end{cases}
```

with the transition blending factor

```math
\gamma = \frac{Re - 2300}{4000 - 2300}.
```

Notes on the modeling choices:

- **Laminar value ``Nu = 4``** is taken as the average of the two classical fully-developed
  constant-property limits, ``4.364`` (uniform heat flux) and ``3.657`` (uniform wall
  temperature) — Eq. 2.42 of Lamarche (2023).
- **Transition region** linearly blends the laminar value with the turbulent correlation
  evaluated at ``Re = 4000``, avoiding the discontinuity of a hard switch.
- ``Nu_{4000}`` is the turbulent expression evaluated at ``Re = 4000``.

See [`Nusselt`](@ref). A convenience overload takes the flow speed and fluid properties
directly and computes ``Re`` and ``Pr`` internally.

### Annulus correction

For the annular channel of a coaxial exchanger the Nusselt number is modified to account for
the curved inner and outer walls. [`Nusselt_annulus`](@ref) follows Lamarche (2021): with the
radius ratio ``a = r_o / r_b``,

```math
Nu_{\text{lam}} = 3.66 + 1.2\sqrt{a},
```

and the turbulent branch applies a geometric factor ``F_a`` and a modified denominator
constant ``k_1`` to the Gnielinski form. The annulus uses an area-weighted equivalent roughness
``\epsilon = (\epsilon_o r_b + \epsilon_i r_o)/(r_b + r_o)`` and the hydraulic radius
``r_b - r_o``.

## Convective resistance

Given the Nusselt number, the convective coefficient and the per-pipe convective resistance are

```math
h = \frac{Nu\,k_f}{2r}, \qquad
R_f = \frac{1}{2\pi r\, h} = \frac{1}{\pi\,Nu\,k_f}.
```

These are Eqs. 2.32 and 5.6 of Lamarche (2023). [`convection_coefficient`](@ref) computes ``h``
alone (either from a pre-computed ``Nu`` or the raw flow speed and properties) — this is what
[`resistance_coaxial`](@ref) uses for its center-pipe and annulus convection coefficients. The
resistance is **per pipe**; the multipole formulas account for the number of pipes in the
borehole. [`resistance_fluid`](@ref) accepts either a pre-computed ``Nu`` or the raw flow speed
and properties.

## Functions on this page

```@docs
Reynolds
Prandtl
friction_factor_Colebrook_White
friction_factor_Tkachenko_Mileikovskyi
Nusselt
Nusselt_annulus
convection_coefficient
resistance_fluid
```
