# API reference

This page is an alphabetical index of every exported symbol. Each entry links to its full
docstring, which lives next to the relevant theory on the modeling pages.

```@index
Modules = [BoreholeResistance]
```

## By topic

### Dimensionless numbers and friction factors

- [`Reynolds`](@ref), [`Prandtl`](@ref)
- [`Nusselt`](@ref), [`Nusselt_annulus`](@ref)
- [`friction_factor_Colebrook_White`](@ref), [`friction_factor_Tkachenko_Mileikovskyi`](@ref)

### Thermal resistances

- [`resistance_fluid`](@ref) — fluid convection
- [`resistance_pipe`](@ref) — pipe wall conduction
- [`resistance_borehole_multipole`](@ref) — borehole resistance ``R_b``
- [`resistance_total_internal_multipole`](@ref) — total internal resistance ``R_a``
- [`resistance_borehole_effective`](@ref) — effective resistance ``R_b^*``

### Water property functions

- [`water_k`](@ref), [`water_cp`](@ref), [`water_ρ`](@ref), [`water_μ`](@ref)

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Fluid convective resistance](@ref), [Pipe conductive resistance](@ref),
    [Borehole (grout) resistance](@ref), [Effective resistance](@ref) and
    [Water properties](@ref).
