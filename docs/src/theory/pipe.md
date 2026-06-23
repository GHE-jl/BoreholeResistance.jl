# Pipe conductive resistance

The pipe wall adds a conductive resistance between the fluid boundary layer and the grout. For
radial conduction through a cylindrical shell of inner radius ``r_i`` and outer radius ``r_o``
made of a material with conductivity ``k_p``, the steady-state resistance per unit length is

```math
R_p = \frac{\ln(r_o / r_i)}{2\pi\,k_p}.
```

This is the textbook cylindrical-shell result (Bergman & Incropera, 2011; Lamarche, 2023). It
is **per pipe** and **independent of the flow rate** — only the geometry and the pipe material
matter. For polyethylene (HDPE) pipe ``k_p \approx 0.4`` W/m·K.

The pipe resistance combines in series with the fluid resistance to form the per-pipe input to
the grout calculation:

```math
R_p^{\text{tot}} = R_f + R_p,
```

as described in the [resistance network](@ref Resistance-network) overview.

## Function on this page

```@docs
resistance_pipe
```
