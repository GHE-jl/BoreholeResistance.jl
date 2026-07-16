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
inside [`Nusselt`](@ref), see [Fluid convective resistance](@ref) for the regime boundaries.

## 4. The individual resistances

```julia
Rp = resistance_pipe(ro, ri, kp)                 # pipe wall conduction
Rf = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ) # fluid convection
```

`Rp` is purely geometric (independent of flow), `Rf` falls as the flow rate, and hence ``Nu``, rises.

## 5. Borehole resistance ``R_b`` and internal resistance ``R_a``

The short form takes the pre-computed `Rp` and `Rf`. `order = 1` selects the first-order
multipole (recommended); `order = 0` is the line-source approximation.

```julia
Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; order = 1)
Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf; order = 1)
```

Equivalently, the **long form** computes `Rf` and `Rp` internally from geometry and flow:

```julia
Rb = resistance_ULoop_borehole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ; order = 1)
```

The grout-only contribution is recovered as ``R_g = R_b - (R_p + R_f)/N``, where ``N`` is the
number of pipes (``N = 2`` for a single U-tube).

## 6. Effective borehole resistance ``R_b^*``

`Rb*` corrects `Rb` for the axial thermal short-circuit between the down-flowing and
up-flowing legs. It needs the flow rate, the borehole length and the volumetric heat capacity:

```julia
Rbe = resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra)
```

or, from geometry directly, the all-in-one form used in the [Quick start](@ref):

```julia
Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf; nLoop = 1)
```

`Rbe ≥ Rb` always — short-circuiting can only degrade performance.

## 7. Double U-tube

Set `nLoop = 2` for a double U-tube (four pipes on a circle of radius). For the
internal resistance you also choose how the two loops are paired with the `network` keyword
(`"diagonal"`, the default, or `"adjacent"`):

```julia
Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop = 2, order = 1)
Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf;
                                         nLoop = 2, order = 1, network = "diagonal")
```

The effective resistance also takes `nLoop = 2`, using the double-U network of Claesson & Javed
(2019). Pass the flow **in one loop** (the total system flow is `2V`) and the `network` that
matches `Ra`:

```julia
Rbe = resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra; nLoop = 2, network = "diagonal")
```

## 8. Coaxial (concentric-tube) exchanger

A coaxial borehole is described by two resistances instead of the U-tube network:
[`resistance_coaxial`](@ref) returns ``R_1`` (annulus → borehole wall, equal to ``R_b``) and
``R_{12}`` (centre pipe ↔ annulus). Pass the four pipe radii and the borehole radius
(``r_b > r_{oo} > r_{oi} > r_{io} > r_{ii}``):

```julia
rii = 0.0247               # inner-pipe inner radius          [m]
rio = 0.0301               # inner-pipe outer radius          [m]
roi = 0.0487               # outer-pipe inner radius          [m]
roo = 0.0572               # outer-pipe outer radius          [m]
kpi = kp                   # inner-pipe conductivity (reuse the HDPE value from §2) [W/m·K]
kpo = kp                   # outer-pipe conductivity                               [W/m·K]

R1, R12 = resistance_coaxial(V, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ)
Rbe = resistance_coaxial_effective(V, H, cf, ρf, R1, R12; model = "UHF")
```

Internally the centre-pipe convection uses [`Nusselt`](@ref) and the annulus uses the
annulus-specific [`Nusselt_annulus`](@ref) (hydraulic radius ``r = r_{oi} - r_{io}``). See
[Borehole (grout) resistance](@ref) for the model and the `model` keyword options.

## Validation scripts

The `script/` directory contains runnable, assertion-checked validation scripts that exercise
every overload across a flow-rate sweep. Run them from the package root:

```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=script/ script/script_single_Uloop.jl
```

| Script | What it validates |
|---|---|
| `script_single_Uloop.jl` | Single U-tube ``R_b``/``R_a`` vs Javed & Spitler (2017), Table 6; plus an HDPE flow sweep. |
| `script_double_Uloop.jl` | Double U-tube ``R_b``, ``R_a`` (diagonal/adjacent) and ``R_b^*`` vs Claesson & Javed (2019), Tables 1–3. |
| `script_fluid_annulus.jl` | Friction factors (Colebrook–White vs Tkachenko–Mileikovskyi), Nusselt pipe vs annulus. |
| `script_coaxial.jl` | Coaxial GHE ``R_1``, ``R_{12}``, ``R_b^*`` (UHF/UBW/mean/gradient) vs Lamarche (2021). |
