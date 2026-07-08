# BoreholeResistance.jl

*Thermal resistances of grouted borehole heat exchangers, in pure Julia.*

`BoreholeResistance.jl` computes the steady-state thermal resistances that govern heat
transfer between the circulating fluid and the ground in a vertical
**ground heat exchanger (GHE)**. It assembles the full fluid-to-ground resistance network:

- the **fluid convective resistance** ``R_f`` (Gnielinski correlation, pipe and annulus),
- the **pipe wall conductive resistance** ``R_p``,
- the **borehole / grout resistance** ``R_b`` (zeroth- and first-order **multipole method**),
- the **total internal resistance** ``R_a`` between the upward and downward legs,
- the **effective borehole resistance** ``R_b^*`` that accounts for axial thermal
  short-circuiting along the borehole depth,

together with temperature-dependent thermophysical properties of water.

The package has **no external dependencies** — everything is implemented with the Julia
standard library — and supports **single U-tube**, **double U-tube** and **coaxial
(concentric-tube)** configurations.

## Why borehole resistance matters

For a borehole carrying a heat load ``q`` [W/m], the temperature difference between the mean
fluid and the borehole wall is

```math
\bar{T}_f - T_b = q \, R_b^* .
```

``R_b^*`` therefore enters directly into every GHE sizing and simulation calculation: a lower
borehole resistance means the fluid runs closer to the ground temperature, improving heat-pump
performance. This package provides the resistances; downstream packages
(see [Ecosystem](@ref)) use them to predict fluid temperatures over time.

## Installation

The package is not yet registered. Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GeothermalJL/BoreholeResistance.jl")
```

or, in the Pkg REPL mode (press `]`):

```
pkg> add https://github.com/GeothermalJL/BoreholeResistance.jl
```

## Quick start

```julia
using BoreholeResistance

# Fluid properties of water at 10 °C
T0 = 10.0
kf = water_k(T0);  cf = water_cp(T0);  ρf = water_ρ(T0);  μf = water_μ(T0)

# Geometry [m] and material conductivities [W/m·K]
H, s, rb, ro, ri = 150.0, 0.05, 0.08, 0.022, 0.017
ks, kg, kp       = 3.0, 1.6, 0.4
V = 30 / 6e4    # 30 L/min expressed in m³/s

# Effective borehole resistance [m·K/W]
Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
```

## Manual outline

- **[Tutorial](@ref)** — a worked single- and double-U-tube example, step by step.
- **Modeling theory** — the physics behind each resistance, with the governing equations and
  their source references:
  - [Resistance network](@ref) — how the individual resistances combine.
  - [Fluid convective resistance](@ref) — ``Re``, ``Pr``, ``Nu``, friction factors, ``R_f``.
  - [Pipe conductive resistance](@ref) — ``R_p``.
  - [Borehole (grout) resistance](@ref) — the multipole method for ``R_b`` and ``R_a``.
  - [Effective resistance](@ref) — ``R_b^*`` and thermal short-circuiting.
- **[Water properties](@ref)** — the polynomial correlations and their validity range.
- **[API reference](@ref)** — the complete docstring reference for every exported function.
- **[References](@ref)** — the bibliography underpinning the implementation.

## Conventions used throughout

| Symbol | Meaning | Unit |
|---|---|---|
| ``R_f`` | Fluid convective resistance (per pipe) | m·K/W |
| ``R_p`` | Pipe wall conductive resistance (per pipe) | m·K/W |
| ``R_b`` | Borehole thermal resistance (fluid → borehole wall) | m·K/W |
| ``R_a`` | Total internal resistance (leg → leg) | m·K/W |
| ``R_b^*`` | Effective borehole resistance (short-circuit corrected) | m·K/W |
| ``k_s, k_g, k_p, k_f`` | Conductivity of ground, grout, pipe, fluid | W/m·K |
| ``r_b, r_o, r_i`` | Borehole, pipe-outer, pipe-inner radius | m |
| ``s`` | Shank spacing (centre-to-centre of the two legs) | m |
| ``H`` | Borehole (active) length | m |
| ``\dot V`` | Mean fluid speed in a pipe | m/s |
| ``V`` | Volumetric flow rate in a pipe | m³/s |

!!! note "Heat-capacity convention"
    `water_cp(T)` returns the **mass-specific** heat ``c_f`` [J/kg·K]. The resistance functions
    take ``c_f`` and ``\rho_f`` separately. The **volumetric** specific heat
    ``C_f = c_f\,\rho_f`` [J/m³·K] is what downstream moving-source models in the ecosystem
    consume — convert explicitly when crossing that boundary.

## Ecosystem

`BoreholeResistance.jl` is the resistance layer of a three-package geothermal stack:

| Package | Role |
|---|---|
| **BoreholeResistance.jl** | Water properties + borehole thermal resistances (this package). |
| **GroundResponse.jl** | Ground *g*-function models and spatial superposition. |
| **GroundHeatExchanger.jl** | Simulation orchestration; depends on and re-exports both. |

Calling `using GroundHeatExchanger` re-exports every function documented here, so the
resistances can be used without importing this package explicitly.
