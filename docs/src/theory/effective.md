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

## Inputs

``R_b^*`` is built from the local resistances ``R_b`` and ``R_a``, the borehole length ``H``,
and the fluid heat-capacity flow rate. The governing quantity is the thermal-capacity resistance
factor

```math
R_V = \frac{H}{V\,c_f\,\rho_f},
```

where ``V`` is the volumetric flow rate **in one U-tube loop** and ``c_f \rho_f`` is the
volumetric heat capacity of the fluid. For a double U-tube the total system flow is ``2V`` (the
two loops in parallel), but ``R_V`` uses the per-loop value ``V`` — pass the per-loop flow, not
the total.

## Two boundary conditions, averaged

The package computes ``R_b^*`` under two idealized boundary conditions — uniform heat flux (UHF)
and uniform borehole-wall temperature (UBW) — and, by default, averages them. The `model`
keyword selects `"UHF"`, `"UBW"` or `"mean"` (the default).

### Single U-tube (`nLoop = 1`)

Following Hellström (1991) / Javed & Spitler (2016), simplified for symmetric legs (Claesson &
Javed, 2019, Eqs. 37–38):

```math
R_{b,\text{UHF}}^* = R_b + \frac{R_V^2}{3 R_a},
\qquad
R_{b,\text{UBW}}^* = R_b\,\eta\coth\eta,
\qquad
\eta = \frac{R_V}{\sqrt{R_b R_a}}.
```

### Double U-tube (`nLoop = 2`)

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

## Overloads

[`resistance_ULoop_effective`](@ref) is available in three forms of increasing convenience,
each accepting `nLoop`, `model` and (for the double U-tube) `network`:

1. from pre-computed ``R_b`` and ``R_a``;
2. from ``R_p`` and ``R_f`` (computes ``R_b`` and ``R_a`` internally);
3. from raw geometry and fluid properties (computes everything).

Both configurations are validated: the single U-tube against Javed & Spitler (2016/2017) and the
double U-tube against Claesson & Javed (2019, Tables 1–3), reproduced by `script_double_Uloop.jl`.

## Function on this page

```@docs
resistance_ULoop_effective
```
