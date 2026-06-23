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
and the fluid heat-capacity flow rate. The governing dimensionless group is the ratio of the
axial advective capacity to the borehole resistance,

```math
\tau = \frac{H}{V\,c_f\,\rho_f},
```

where ``V`` is the volumetric flow rate per pipe and ``c_f \rho_f`` is the volumetric heat
capacity of the fluid.

## Delta-network mapping

The two legs and the borehole wall form a delta network (see
[Resistance network](@ref Resistance-network)). The leg-to-wall and leg-to-leg resistances are

```math
R_1 = 2 R_b, \qquad
R_{12} = \frac{2 R_a R_1}{2 R_1 - R_a}.
```

## Two boundary conditions, averaged

The package follows Javed & Spitler (2016) and computes ``R_b^*`` under two idealized boundary
conditions, then averages them.

### Uniform borehole wall temperature (UBW)

Define

```math
\eta = \frac{\tau}{2 R_b}\sqrt{1 + \frac{4 R_b}{R_{12}}}.
```

Then

```math
R_{b,\text{UBW}}^* =
\begin{cases}
R_b\,\eta\,\coth\eta & \eta > 1 \\[6pt]
R_b + \dfrac{\tau^2}{3 R_{12}} + \dfrac{\tau^2}{12 R_b} & \eta \le 1
\end{cases}
```

(Eqs. 3.68–3.70 of Javed & Spitler, 2016). The small-``\eta`` branch is the series expansion
that stays numerically well-behaved when the short-circuit is weak.

### Uniform heat flux (UHF)

```math
R_{b,\text{UHF}}^* = R_b + \frac{\tau^2}{3 R_a}
```

(Eq. 3.67 of Javed & Spitler, 2016).

### Average

The reported effective resistance is the mean of the two limits, which is the recommended
practical estimate:

```math
R_b^* = \tfrac{1}{2}\left(R_{b,\text{UBW}}^* + R_{b,\text{UHF}}^*\right).
```

## Overloads

[`resistance_borehole_effective`](@ref) is available in three forms of increasing convenience:

1. from pre-computed ``R_b`` and ``R_a``;
2. from ``R_p`` and ``R_f`` (computes ``R_b`` and ``R_a`` internally);
3. from raw geometry and fluid properties (computes everything).

!!! warning "Single U-tube derivation"
    The effective-resistance formulas are derived for the single U-tube (`nLoop = 1`, two pipes
    per borehole). Applying them to a double U-tube is an approximation.

## Function on this page

```@docs
resistance_borehole_effective
```
