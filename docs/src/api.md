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

- [`convection_coefficient`](@ref) — convective heat transfer coefficient ``h``
- [`resistance_fluid`](@ref) — fluid convection ``R_f``
- [`resistance_pipe`](@ref) — pipe wall conduction ``R_p``
- [`resistance_ULoop_borehole`](@ref) — borehole resistance ``R_b``
- [`resistance_ULoop_total_internal`](@ref) — total internal resistance ``R_a``
- [`resistance_ULoop_effective`](@ref) — effective resistance ``R_b^*``
- [`resistance_coaxial`](@ref) — coaxial resistances ``R_1``, ``R_{12}``
- [`resistance_coaxial_effective`](@ref) — coaxial effective resistance ``R_b^*``

### Water property functions

- [`water_k`](@ref), [`water_cp`](@ref), [`water_ρ`](@ref), [`water_μ`](@ref)
- [`fluid_property`](@ref) — water/MPG/MEG properties via CoolProp

!!! tip "Where the docstrings live"
    Full signatures and argument lists are rendered inline on the theory pages:
    [Fluid convective resistance](@ref), [Pipe conductive resistance](@ref),
    [Borehole (grout) resistance](@ref), [Effective resistance](@ref) and
    [Water properties](@ref).
