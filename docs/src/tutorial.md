# Tutorial

This tutorial walks through a complete borehole-resistance calculation, from fluid properties
to the effective borehole resistance ``R_b^*``, for both a single and a double U-tube. Every
function used here is documented in the [API reference](@ref).

## 1. Fluid properties

All correlations take the mean fluid temperature in degrees Celsius and return SI units. Pick
a representative loop temperature (here 10 °C) and evaluate the four properties needed by the
resistance functions:

```julia
using BoreholeResistance

T0 = 10.0                 # mean fluid temperature [°C]
kf = water_k(T0)          # thermal conductivity   [W/m·K]
cf = water_cp(T0)         # specific heat           [J/kg·K]
ρf = water_ρ(T0)          # density                 [kg/m³]
μf = water_μ(T0)          # dynamic viscosity       [Pa·s]
```

See [Water properties](@ref) for the underlying correlations and their validity range.

## 2. Geometry and flow

```julia
H  = 100.0                # active borehole length          [m]
rb = 0.075                # borehole radius                 [m]
ro = 0.020                # pipe outer radius               [m]
ri = 0.0164               # pipe inner radius               [m]
s  = 0.08                 # shank spacing (centre-to-centre)[m]

ks = 2.0                  # ground conductivity             [W/m·K]
kg = 1.0                  # grout conductivity              [W/m·K]
kp = 0.4                  # pipe (HDPE) conductivity        [W/m·K]
ϵ  = 5e-6                 # pipe roughness                  [m]

V  = 15.0 / 1000 / 60     # volumetric flow per pipe        [m³/s]  (15 L/min)
V̇ = V / (π * ri^2)       # mean fluid speed                [m/s]
```

## 3. Flow regime: Reynolds, Prandtl, Nusselt

```julia
Re = Reynolds(V̇, ri, ρf, μf)
Pr = Prandtl(kf, cf, μf)
Nu = Nusselt(Re, Pr, ri, ϵ)        # Gnielinski correlation
```

The transition through laminar / transitional / turbulent regimes is handled automatically
inside [`Nusselt`](@ref); see [Fluid convective resistance](@ref) for the regime boundaries.

## 4. The individual resistances

```julia
Rp = resistance_pipe(ro, ri, kp)                 # pipe wall conduction
Rf = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ) # fluid convection
```

`Rp` is purely geometric (independent of flow); `Rf` falls as the flow rate — and hence ``Nu``
— rises.

## 5. Borehole resistance ``R_b`` and internal resistance ``R_a``

The short form takes the pre-computed `Rp` and `Rf`. `order = 1` selects the first-order
multipole (recommended); `order = 0` is the line-source approximation.

```julia
Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf; order = 1)
Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf; order = 1)
```

Equivalently, the **long form** computes `Rf` and `Rp` internally from geometry and flow:

```julia
Rb = resistance_borehole_multipole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ; order = 1)
```

The grout-only contribution is recovered as ``R_g = R_b - R_p - R_f``.

## 6. Effective borehole resistance ``R_b^*``

`Rb*` corrects `Rb` for the axial thermal short-circuit between the down-flowing and
up-flowing legs. It needs the flow rate, the borehole length and the volumetric heat capacity:

```julia
Rbe = resistance_borehole_effective(V, H, cf, ρf, Rb, Ra)
```

or, from geometry directly, the all-in-one form used in the [Quick start](@ref BoreholeResistance.jl):

```julia
Rbe = resistance_borehole_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf; nLoop = 1)
```

`Rbe ≥ Rb` always — short-circuiting can only degrade performance.

## 7. Double U-tube

Set `nLoop = 2` for a double U-tube (four pipes). For the internal resistance you can also
choose how the legs are paired with the `network` keyword (`"diagonal"`, the default, or
`"adjacent"`):

```julia
Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf; nLoop = 2, order = 1)
Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf;
                                         nLoop = 2, order = 1, network = "diagonal")
```

!!! warning "Effective resistance is single-U only"
    `resistance_borehole_effective` is derived for the single U-tube (`nLoop = 1`).
    Use it with double-U resistances only as an approximation.

## 8. Annulus / coaxial helpers

For the annular channel of a coaxial exchanger, the Reynolds number and friction factors take
the hydraulic radius ``r = r_b - r_o``, and the Nusselt number uses the annulus-specific
correlation:

```julia
Nu_a = Nusselt_annulus(V̇, rb, ro, kf, cf, ρf, μf)
```

The full coaxial borehole resistance is not yet implemented (see the project TODO).

## Validation scripts

The `script/` directory contains runnable, assertion-checked validation scripts that exercise
every overload across a flow-rate sweep. Run them from the package root:

```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=script/ script/script_single_Uloop.jl
```

| Script | What it validates |
|---|---|
| `script_single_Uloop.jl` | Single U-tube: `Re`, `Pr`, `Nu`, `Rf`, `Rp`, `Rb` (order 0/1), `Ra`, `Rb*`. |
| `script_double_Uloop.jl` | Double U-tube: `Rb`, `Ra` diagonal/adjacent, `Rbe`. |
| `script_annulus.jl` | Friction factors CW vs TM, Nusselt pipe vs annulus. |
| `script_coaxial.jl` | Coaxial GHE (placeholder — not yet implemented). |
