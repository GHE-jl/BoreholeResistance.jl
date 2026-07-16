# Resistance network

This page describes the physical picture behind the package: how heat travels from the
circulating fluid to the surrounding ground, and how that path is decomposed into a network of
thermal resistances. The subsequent theory pages derive each resistance in turn.

## The borehole heat exchanger

A vertical GHE is a borehole of radius ``r_b`` drilled to a depth ``H``, into which one or more
U-shaped pipe loops or a coaxial configuration are inserted and the remaining space backfilled with **grout**. A heat
carrier fluid (water, or water with antifreeze) is pumped down one leg of the U-tube and back
up the other, exchanging heat with the ground through the pipe wall and the grout.

The radial heat path, from the inside out, crosses four media:

1. the **fluid** boundary layer (convection),
2. the **pipe wall** (conduction),
3. the **grout** filling the borehole (conduction, multi-pipe geometry),
4. the **ground** outside the borehole wall (handled separately by transient *g*-function models).

This package computes the steady-state resistance from the fluid up to the **borehole wall**
at ``r = r_b``. Everything beyond the wall, the transient response of the ground, is the domain of the `GroundResponse.jl` package.

## Definition of the borehole resistance

The **borehole thermal resistance** ``R_b`` [m·K/W] relates the heat transfer rate per unit
borehole length ``q`` [W/m] to the temperature difference between the mean fluid and the
borehole wall:

```math
\bar{T}_f - T_b = q R_b .
```

It is a *local*, steady-state quantity: it assumes the heat flux is uniform along the depth and
ignores how the fluid temperature itself changes as it travels down and up the borehole. That
last effect, the **thermal short-circuit** between the two legs, is captured separately by
the effective resistance ``R_b^*``.

## Series decomposition

For a single pipe, the resistance from the fluid to the outer pipe wall is a simple series sum
of the convective and conductive parts:

```math
R_p^{\text{tot}} = R_f + R_p,
```

where

- ``R_f`` is the **fluid convective resistance** of one pipe — see
  [Fluid convective resistance](@ref);
- ``R_p`` is the **pipe wall conductive resistance** of one pipe — see
  [Pipe conductive resistance](@ref).

This combined per-pipe resistance ``R_p^{\text{tot}}`` is the input to the grout calculation.
Throughout the multipole formulas it appears through the dimensionless group

```math
\beta = 2\pi k_g \, R_p^{\text{tot}} = 2\pi k_g (R_f + R_p),
```

which measures the pipe-plus-fluid resistance relative to the grout conductivity ``k_g``.

## The grout step is not a simple series resistance

The grout contribution cannot be written as a plain series resistance, because the pipes sit
**off-centre** inside the borehole and there are **several of them** (two for a single,
four for a double U-tube). The temperature field in the grout is genuinely two-dimensional, and the
legs thermally interact.

The **multipole method** of Hellström (1991) solves this 2-D conduction problem by expanding
the grout temperature field in a series of multipoles (line source + dipole + quadrupole + …).
Truncating the expansion at:

- **order 0** recovers the classical line-source result;
- **order 1** adds the first multipole correction and is accurate to better than ~1 % for
  typical geometries.

Both orders are implemented; see [Borehole (grout) resistance](@ref). The method also yields
the **total internal resistance** ``R_a`` between the up and down legs, which is what links the
local ``R_b`` to the depth-corrected ``R_b^*``.

## The delta network

The two legs of a single U-tube and the borehole wall form a three-node **delta (Δ) network**
of three resistances:

```math
R_1 = R_2 \quad (\text{leg-to-wall, equal by symmetry}), \qquad R_{12}\ (\text{leg-to-leg}).
```

The multipole outputs map onto this network as

```math
R_1 = 2 R_b, \qquad R_{12} = \frac{2 R_a R_1}{2 R_1 - R_a}.
```

This network underlies the [effective resistance](@ref Effective-resistance) calculation, which
folds in the axial fluid-temperature variation along ``H``. For the symmetric U-tube it reduces
to the compact closed forms in ``R_b`` and ``R_a`` given on that page (the double U-tube uses the
analogous four-pipe network of Claesson & Javed, 2019).

The **coaxial** configuration does not fit this delta picture: there is no off-centre multipole
problem, since the two flow channels (centre pipe and annulus) are concentric. Instead it reduces
to a simple two-resistance series network, ``R_{12}`` (centre pipe ↔ annulus) and ``R_1`` (annulus
→ borehole wall), described in [Coaxial (concentric-tube) exchanger](@ref) and its own effective
resistance in [Effective coaxial resistance](@ref).

## Map of the implementation

| Quantity | Symbol | Function | Theory page |
|---|---|---|---|
| Reynolds, Prandtl, Nusselt | ``Re,Pr,Nu`` | [`Reynolds`](@ref), [`Prandtl`](@ref), [`Nusselt`](@ref) | [Fluid](@ref Fluid-convective-resistance) |
| Friction factor | ``f`` | [`friction_factor_Colebrook_White`](@ref) | [Fluid](@ref Fluid-convective-resistance) |
| Fluid convective resistance | ``R_f`` | [`resistance_fluid`](@ref) | [Fluid](@ref Fluid-convective-resistance) |
| Pipe conductive resistance | ``R_p`` | [`resistance_pipe`](@ref) | [Pipe](@ref Pipe-conductive-resistance) |
| Borehole resistance | ``R_b`` | [`resistance_ULoop_borehole`](@ref) | [Borehole](@ref Borehole-(grout)-resistance) |
| Total internal resistance | ``R_a`` | [`resistance_ULoop_total_internal`](@ref) | [Borehole](@ref Borehole-(grout)-resistance) |
| Effective resistance | ``R_b^*`` | [`resistance_ULoop_effective`](@ref) | [Effective](@ref Effective-resistance) |
| Coaxial resistances | ``R_1,R_{12}`` | [`resistance_coaxial`](@ref) | [Borehole](@ref Borehole-(grout)-resistance) |
| Coaxial effective resistance | ``R_b^*`` | [`resistance_coaxial_effective`](@ref) | [Effective](@ref Effective-resistance) |
