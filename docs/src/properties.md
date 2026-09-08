# Fluid properties

The resistance correlations need the thermophysical properties of the heat carrier fluid as a
function of temperature. [`fluid_property`](@ref) is the recommended entry point: given a
temperature in degrees Celsius, a fluid symbol, and (for mixtures) a mass concentration, it
returns all four properties at once via [CoolProp](https://coolprop.org/fluid_properties/Incompressibles.html):

```julia
using BoreholeResistance

T0 = 10.0
kf, cf, ρf, μf = fluid_property(T0, :water)   # k ≈ 0.578 W/m·K, cp ≈ 4192 J/kg·K, ρ ≈ 999.7 kg/m³, μ ≈ 1.3e-3 Pa·s
```

The `fluid` symbol is `:water` for pure water, or one of `:MPG`, `:MEG`, `:MMA`, `:MEA`, `:MKA`,
`:MKF` for an aqueous mixture (propylene glycol, ethylene glycol, methanol, ethanol, potassium
acetate, or potassium formate respectively), with `percentage` giving the mass fraction of the
additive [%m]:

```julia
k, cp, ρ, μ = fluid_property(0.0, :MPG; percentage = 30)
```

Temperature and concentration validity are fluid- and composition-dependent (e.g. a glycol
mixture's freezing point depends on its concentration); CoolProp raises an error for an
out-of-range mixture request rather than returning a silently wrong value. Pure water has no such
built-in guard, so `fluid_property` checks against the actual boiling point at the reference
pressure and warns for `:water` outside the liquid range.

## Heat-capacity convention

`fluid_property` returns the **mass-specific** heat ``c_f`` [J/kg·K] as its second output. The
resistance functions in this package take ``c_f`` and ``\rho_f`` as separate arguments. Models
that need the **volumetric** specific heat ``C_f`` [J/m³·K], notably the moving-source models in
`GroundResponse.jl` and the fluid-temperature routines in `GroundHeatExchanger.jl`, require the
explicit product:

```math
C_f = c_f \rho_f.
```

Keep track of which convention a downstream function expects when crossing package boundaries.

## Legacy polynomial fits (deprecated)

Before [`fluid_property`](@ref) was added, the package shipped standalone polynomial fits for
**pure water only**, valid from 0 °C to 100 °C at near-atmospheric pressure, fitted to data from
the Engineering ToolBox:

| Function | Property | Symbol | Unit |
|---|---|---|---|
| [`water_k`](@ref) | Thermal conductivity | ``k_f`` | W/m·K |
| [`water_cp`](@ref) | Specific heat capacity | ``c_f`` | J/kg·K |
| [`water_ρ`](@ref) | Density | ``\rho_f`` | kg/m³ |
| [`water_μ`](@ref) | Dynamic viscosity | ``\mu_f`` | Pa·s |

These four functions are **deprecated** in favor of `fluid_property(T, :water)` and emit a
deprecation warning on every call. They are kept, still exported, and numerically unchanged for
backward compatibility (they are also re-exported by `GroundHeatExchanger.jl`), but new code
should use [`fluid_property`](@ref) instead.

## Functions on this page

```@docs
fluid_property
water_k
water_cp
water_ρ
water_μ
```
