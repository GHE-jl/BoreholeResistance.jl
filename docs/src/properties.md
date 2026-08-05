# Water properties

The resistance correlations need the thermophysical properties of the heat carrier fluid as a
function of temperature. The package ships polynomial fits for **pure water**, valid from 0 °C
to 100 °C at near-atmospheric pressure, fitted to data from the Engineering ToolBox. Each
function takes the temperature in degrees Celsius and returns SI units, and emits a warning if
called outside the fitted range.

| Function | Property | Symbol | Unit |
|---|---|---|---|
| [`water_k`](@ref) | Thermal conductivity | ``k_f`` | W/m·K |
| [`water_cp`](@ref) | Specific heat capacity | ``c_f`` | J/kg·K |
| [`water_ρ`](@ref) | Density | ``\rho_f`` | kg/m³ |
| [`water_μ`](@ref) | Dynamic viscosity | ``\mu_f`` | Pa·s |

Each fit is a degree-5 polynomial,

```math
\phi(T) = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4 + a_5 T^5,
```

reproducing the reference data to better than ~0.5 % over the validity range (1 % for
viscosity).

## Heat-capacity convention

`water_cp` returns the **mass-specific** heat ``c_f`` [J/kg·K]. The resistance functions in this
package take ``c_f`` and ``\rho_f`` as separate arguments. Models that need the **volumetric**
specific heat ``C_f`` [J/m³·K], notably the moving-source models in `GroundResponse.jl` and
the fluid-temperature routines in `GroundHeatExchanger.jl`, require the explicit product:

```math
C_f = c_f \rho_f = \texttt{water\_cp(T)} \times \texttt{water\_ρ(T)}.
```

Keep track of which convention a downstream function expects when crossing package boundaries.

## Example

```julia
using BoreholeResistance

T0 = 10.0
kf = water_k(T0)          # ≈ 0.578  W/m·K
cf = water_cp(T0)         # ≈ 4192   J/kg·K
ρf = water_ρ(T0)          # ≈ 999.7  kg/m³
μf = water_μ(T0)          # ≈ 1.3e-3 Pa·s
Cf = cf * ρf              # volumetric specific heat [J/m³·K]
```

## Antifreeze mixtures (CoolProp)

The polynomial fits above only cover pure water. For antifreeze mixtures — aqueous propylene
glycol (`:MPG`) or ethylene glycol (`:MEG`) — [`fluid_property`](@ref) wraps
[CoolProp](https://coolprop.org/fluid_properties/Incompressibles.html) to return all four
properties (``k``, volumetric heat capacity, ``\rho``, ``\mu``) at once, given a temperature in
kelvin and an optional mass concentration.

```julia
using BoreholeResistance

k, Cs, ρ, μ = fluid_property(273.15, :MPG; percentage = 30)
```

## Functions on this page

```@docs
water_k
water_cp
water_ρ
water_μ
fluid_property
```
