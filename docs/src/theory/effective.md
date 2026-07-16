# Effective resistance

The local borehole resistance ``R_b`` assumes a uniform heat flux along the depth. In reality
the fluid changes temperature as it descends one leg and rises in the other, so the down-leg and
up-leg exchange heat with each other through the grout — the **thermal short-circuit**. The
**effective borehole resistance** ``R_b^*`` (also written ``R_b^{\text{eff}}`` or ``R_{be}``)
folds this axial effect into a single depth-averaged resistance, and is the value that should
be used in sizing and simulation.

Because the short-circuit always adds a parasitic heat path,

```math
R_b^* \ge R_b .
```

The correction grows with borehole length ``H`` and shrinks with flow rate — faster flow means
the fluid spends less time exchanging heat with the opposite leg.

``R_b^*`` is built from the local resistances ``R_b`` and ``R_a``, the borehole length ``H``,
and the fluid heat-capacity flow rate. The governing quantity is the thermal-capacity resistance
factor

```math
R_V = \frac{H}{VC_f},
```

where ``V`` is the volumetric flow rate **in one U-tube loop** and ``C_f`` is the
volumetric heat capacity of the fluid. For a double U-tube the total system flow is ``2V`` (the
two loops in parallel), but ``R_V`` uses the per-loop value ``V`` — pass the per-loop flow, not
the total.

## Effective U-loop resistance

### Two boundary conditions, averaged

The package computes ``R_b^*`` under two idealized boundary conditions — uniform heat flux (UHF)
and uniform borehole-wall temperature (UBW) — and, by default, averages them. The `model`
keyword selects `"UHF"`, `"UBW"` or `"mean"` (the default).

#### Single U-tube (`nLoop = 1`)

Following Hellström (1991) / Javed & Spitler (2016), simplified for symmetric legs (Claesson &
Javed, 2019, Eqs. 37–38):

```math
R_{b,\text{UHF}}^* = R_b + \frac{R_V^2}{3 R_a},
\qquad
R_{b,\text{UBW}}^* = R_b\,\eta\coth\eta,
\qquad
\eta = \frac{R_V}{\sqrt{R_b R_a}}.
```

#### Double U-tube (`nLoop = 2`)

The two loops make the internal coupling stronger, changing the network coefficients (Claesson &
Javed, 2019, Eqs. 44 and 46):

```math
R_{b,\text{UHF}}^* = R_b + \frac{R_V^2}{6 R_a},
\qquad
R_{b,\text{UBW}}^* = R_b\,\eta\coth\eta,
\qquad
\eta = \frac{R_V}{\sqrt{2 R_b R_a}}.
```

For the double U-tube, ``R_a`` must be the internal resistance of the matching flow
configuration — pass the same `network` (`"diagonal"` or `"adjacent"`) that was used for
``R_a``.

### Average

The default `model = "mean"` returns

```math
R_b^* = \tfrac{1}{2}\left(R_{b,\text{UHF}}^* + R_{b,\text{UBW}}^*\right).
```

### Overloads

[`resistance_ULoop_effective`](@ref) is available in three forms of increasing convenience,
each accepting `nLoop`, `model` and (for the double U-tube) `network`:

1. from pre-computed ``R_b`` and ``R_a``;
2. from ``R_p`` and ``R_f`` (computes ``R_b`` and ``R_a`` internally);
3. from raw geometry and fluid properties (computes everything).

Both configurations are validated: the single U-tube against Javed & Spitler (2016/2017) and the
double U-tube against Claesson & Javed (2019, Tables 1–3), reproduced by `script_double_Uloop.jl`.

## Effective coaxial resistance

With the groups (Eq. 7, ``\dot m c_f = VC_f`` of Lamarche, 2021)

```math
\gamma = \frac{H}{2\,\dot VC_f\,R_1}, \qquad
\xi = \sqrt{\frac{R_a}{4 R_1}}, \qquad
\eta = \frac{\gamma}{\xi},
```

the two closed forms exposed by [`resistance_coaxial_effective`](@ref) are

```math
R_b^* = R_1\,\eta\coth\eta \qquad (\text{UBW, Eq. 14}),
\qquad
R_b^* = R_1\left(1 + \frac{R_a}{R_{12}}\frac{\eta^2}{3}\right) \qquad (\text{UHF, Eq. 31}).
```

Both flow directions ("center-in" and "annulus-in") give the same ``R_b^*``. Lamarche (2021)
recommends the uniform-heat-flux form (Eq. 31) as the better compromise for coaxial exchangers,
which is the default `model = "UHF"`.

### Linearly-varying far-field temperature (`"UHF_gradient"`)

The models above assume a uniform far-field temperature. When it instead varies linearly with
depth — a geothermal gradient in a deep borehole — Lamarche (2021) shows (Section 3) that a
linearly-varying borehole-wall temperature is the physically consistent extension, but its closed
form (Eq. 43) needs the *actual* wall-temperature slope and reference temperatures as extra
inputs, breaking the borehole/ground decoupling that keeps `"UHF"`/`"UBW"`/`"mean"`
self-contained. Section 4.1 instead proposes a much cheaper proxy: keep the borehole decoupled
from the ground, but let the heat flux itself vary linearly along the borehole instead of being
uniform (Eq. 58), with boundary conditions chosen so it vanishes at one end:

```math
\tilde q'(\tilde z) = 2(1-\tilde z) \quad \text{(heat injection, annulus-in)}, \qquad
\tilde q'(\tilde z) = 2\tilde z \quad \text{(heat extraction, center-in)}.
```

Both cases integrate (Eqs. 60–63) to the same closed form, exposed as `model = "UHF_gradient"`:

```math
R_b^* = R_1\left(1 + \frac{H}{6\,\dot m c_f R_1} + \frac{H^2}{4\,(\dot m c_f)^2 R_1 R_{12}}\right)
\quad (\text{Eqs. 61/63}).
```

It needs no numeric input beyond `"UHF"` (same `V`, `H`, `R_1`, `R_{12}`), but it is **not** a
generalization of `"UHF"` — it is a fixed linear heat-flux shape, not a tunable gradient
magnitude, so it does not reduce to Eq. 31 when there happens to be no gradient (Table 3 of
Lamarche 2021 reports both values for the same borehole, and they differ: `0.0356` mK/W for
`"UHF"` versus `0.0420`/`0.0494` mK/W for `"UHF_gradient"`). It is also only valid for the
flow-direction/heat-mode pairing that Lamarche (2021) identifies as *unfavorable* for the
gradient's sign — heat injection with "annulus-in" or heat extraction with "center-in" when the
far-field temperature increases with depth (mirror the pairing for a negative gradient). For the
*favorable* pairing, Lamarche (2021) found that plain `"UHF"` remains the better estimate; picking
the right model for the situation is the caller's responsibility, since flow direction, heat mode
and gradient sign are not arguments of [`resistance_coaxial_effective`](@ref).

## Function on this page

```@docs
resistance_ULoop_effective
resistance_coaxial_effective
```
