# BoreholeResistance.jl

A Julia package for computing borehole thermal resistances in ground heat exchanger (GHE) systems. Provides fluid, pipe, and grout resistances using the multipole method, along with water thermophysical property functions.

## Quick start

```julia
using BoreholeResistance

T0 = 10.0
kf = water_k(T0);  cf = water_cp(T0);  ρf = water_ρ(T0);  μf = water_μ(T0)

H, s, rb, ro, ri = 150.0, 0.05, 0.08, 0.022, 0.017
ks, kg, kp = 3.0, 1.6, 0.4
V = 30/6e4   # 30 L/min in m³/s

Rb  = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
```

## Water thermophysical properties

Polynomial fits to Engineering Toolbox data, valid 0–100°C at 1 atm:

| Function | Output | Unit |
|---|---|---|
| `water_k(T)` | Thermal conductivity | W/m·K |
| `water_cp(T)` | Specific heat capacity | J/kg·K |
| `water_ρ(T)` | Density | kg/m³ |
| `water_μ(T)` | Dynamic viscosity | Pa·s |

**Convention:** `water_cp` returns mass-specific heat `cf` [J/kg·K]. The volumetric specific heat
`Cf = water_cp(T) * water_ρ(T)` [J/m³·K] is needed separately for `outlet_temperature` /
`inlet_temperature` in GroundHeatExchanger.jl. Borehole resistance functions take `cf` and `ρf`
individually.

## Dimensionless numbers and friction factors

```julia
Reynolds(V̇, r, ρf, μf)                            # V̇: fluid speed [m/s], r: pipe radius [m]
Prandtl(kf, cf, μf)                                # or Prandtl(kf, Cf, ρf, μf) with Cf [J/m³K]
Nusselt(Re, Pr, r, ϵ=5e-6)                         # Gnielinski (laminar / transition / turbulent)
Nusselt_annulus(Re, Pr, rb, ro, ϵo=5e-6, ϵi=5e-6)  # annulus correction — Lamarche (2021)
friction_factor_Colebrook_White(Re, r, ϵ)           # iterative, Darcy friction factor
friction_factor_Tkachenko_Mileikovskyi(Re, r, ϵ)    # explicit approximation (agrees within 1%)
```

## Thermal resistances

All resistances are in [m·K/W].

### Fluid convective resistance
```julia
resistance_fluid(Nu, r, kf)
resistance_fluid(V̇, r, kf, cf, ρf, μf, ϵ=5e-6)
```

### Pipe wall conductive resistance
```julia
resistance_pipe(ro, ri, kp)    # log(ro/ri) / (2π·kp)
```

### Borehole resistance — multipole method (Hellström 1991; Javed & Spitler 2017; Claesson & Javed 2019)
```julia
# Short form: pass pre-computed Rp and Rf
resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=1, order=1)

# Long form: compute Rp/Rf internally from pipe geometry and flow
resistance_ULoop_borehole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0;
    nLoop=1, order=1)
```

`nLoop=1` → single U-tube (2 pipes); `nLoop=2` → double U-tube (4 pipes on a circle of radius
`rc = s/2`, i.e. `s` is the distance between diagonally opposite pipes).
`order=0` → zeroth-order (line-source); `order=1` → first-order multipole.

### Total internal resistance (for Rb* computation)
```julia
resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf;
    nLoop=1, order=1, network="diagonal")
# network="adjacent" for double U-tube with adjacent pipe pair connections
```

### Effective borehole resistance Rb*
Accounts for axial thermal short-circuit along the borehole depth. `V` is the flow rate **in one
U-tube loop** (for a double U-tube the total system flow is `2V`). `model` selects the boundary
condition (`"UHF"`, `"UBW"` or `"mean"`, default `"mean"`); for the double U-tube, `network`
must match the one used for `Ra`:
```julia
# From pre-computed Rb and Ra
resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra; nLoop=1, model="mean")

# From Rp/Rf
resistance_ULoop_effective(V, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf;
    nLoop=1, model="mean", network="diagonal")

# Full form: all geometry and fluid properties
resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0;
    nLoop=1, model="mean", network="diagonal")
```
Single-U formulas follow Javed & Spitler (2016); double-U (`nLoop=2`) follows Claesson & Javed
(2019, Eqs. 44/46).

### Coaxial (concentric-tube) borehole — Lamarche (2021)

For coaxial exchangers the borehole is described by two resistances instead of the U-tube
multipole network: `R12` between the center pipe and the annulus (Eq. 1), and `R1` between the
annulus fluid and the borehole wall (Eq. 2). By Eq. 8, `R1` is also the (steady) borehole
resistance `Rb`.

```julia
# Two resistances (R1, R12) from convection coefficients (hin in the center pipe, hann in annulus)
R1, R12 = resistance_coaxial(rii, rio, roi, roo, rb, kg, kpi, kpo, hin, hann)

# Or computed internally from flow rate and fluid properties (Nusselt / Nusselt_annulus)
R1, R12 = resistance_coaxial(V, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ=0.0)
```

Radii: `rii`/`rio` inner pipe inner/outer, `roi`/`roo` outer pipe inner/outer, `rb` borehole.

The effective resistance `Rb*` accounts for the axial short-circuit between center and annulus.
The `model` keyword selects the closed form: `"UHF"` (uniform heat flux, Eq. 31 — recommended for
coaxial), `"UBW"` (uniform wall temperature, Eq. 14), `"mean"`, or `"UHF_gradient"` (Eqs. 61/63 —
for a linearly-varying far-field temperature, see caveats below). Both flow directions
(center-in / annulus-in) give the same `Rb*` for `"UHF"`/`"UBW"`/`"mean"`.

```julia
# From pre-computed R1 and R12
resistance_coaxial_effective(V, H, cf, ρf, R1, R12; model="UHF")

# Full form: all geometry and fluid properties
resistance_coaxial_effective(V, H, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf,
    ϵ=0.0; model="UHF")
```

`"UHF_gradient"` needs no extra numeric input, but is only valid for the flow-direction/heat-mode
pairing that is *unfavorable* for the far-field gradient's sign (e.g. heat injection with
"annulus-in" when the far-field temperature increases with depth) — see the
[`resistance_coaxial_effective`](https://github.com/GHE-jl/BoreholeResistance.jl/blob/main/src/resistance_borehole.jl)
docstring for the full applicability rule. It is a distinct closed form, not a generalization of
`"UHF"`, and does not reduce to it in the absence of a gradient.

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`.
First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

| Script | What it validates |
|---|---|
| `script_single_Uloop.jl` | Single U-tube Rb/Ra/Rg vs Javed & Spitler (2017) Table 6; HDPE flow sweep; all overloads |
| `script_double_Uloop.jl` | Double U-tube Rb, Ra (diagonal/adjacent), Rb* vs Claesson & Javed (2019) Tables 1–3; all overloads |
| `script_fluid_annulus.jl` | Friction factors CW vs TM, Nusselt pipe vs annulus, all overloads |
| `script_coaxial.jl` | Coaxial GHE: R1, R12, Rb* (UHF/UBW/mean/gradient) vs Lamarche (2021) Boreholes 1 & 3b |

## Installation

The package is not yet registered. Install directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/BoreholeResistance.jl")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add https://github.com/GHE-jl/BoreholeResistance.jl
```

## Dependencies

No external dependencies. All required functions are implemented in pure Julia using only the
standard library.

## Integration with GroundHeatExchanger.jl

`BoreholeResistance.jl` is a dependency of `GroundHeatExchanger.jl`. All resistance functions
and water property functions are re-exported from the top-level package:

```julia
using GroundHeatExchanger   # pulls in BoreholeResistance + GroundResponse

cf = water_cp(T0)           # specific heat [J/kg·K]
ρf = water_ρ(T0)            # density [kg/m³]
Cf = cf * ρf                # volumetric specific heat [J/m³·K]

Rb = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)

# head_loss_Darcy_Weisbach is defined in GroundHeatExchanger.jl (hydraulic / pump sizing)
```

## References

- Mileikovskyi, V., & Tkachenko, T. (2021). Precise Explicit Approximations of the
  Colebrook-White Equation for Engineering Systems. In Z. Blikharskyy (Ed.), Proceedings of 
  EcoComfort 2020 (pp. 303–310). Springer International Publishing.
  https://doi.org/10.1007/978-3-030-57340-9_37
- Hellström, G. (1991). Ground Heat Storage: Thermal Analyses of Duct Storage Systems.
  Lund University.
- Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
  for grouted single U-tube GHEs. Applied Energy, 187, 790–806.
- Claesson, J., & Javed, S. (2019). Explicit multipole formulas for double U-pipe borehole
  heat exchangers. Science and Technology for the Built Environment, 25(8), 980–992.
- Lamarche, L. (2021). Analytic models and effective resistances for coaxial GHEs.
  Geothermics, 97, 102224.
- Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application.
  Springer Nature Switzerland.
