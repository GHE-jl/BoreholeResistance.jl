# Borehole (grout) resistance

The borehole resistance ``R_b`` combines the fluid, pipe and grout contributions into the
single resistance from the fluid to the borehole wall. Because the pipes sit off-centre and
interact, the grout step is solved with the **multipole method** of Hellström (1991), in the
explicit form given by Javed & Spitler (2017) for single U-tubes and Claesson & Javed (2019)
for double U-tubes. The same machinery also yields the **total internal resistance** ``R_a``
between the two legs. Both are implemented in
[`resistance_borehole.jl`](https://github.com/GeothermalJL/BoreholeResistance.jl/blob/master/src/resistance_borehole.jl).

## Common dimensionless groups

Every multipole formula is built from three groups. With the combined per-pipe resistance
``R_p^{\text{tot}} = R_f + R_p``:

```math
\beta = 2\pi k_g \, R_p^{\text{tot}},
\qquad
\sigma = \frac{k_g - k_s}{k_g + k_s},
```

- ``\beta`` scales the pipe-plus-fluid resistance against the grout conductivity;
- ``\sigma`` is the **grout–ground conductivity contrast**, ranging from ``-1`` (highly
  conductive grout) to ``+1`` (insulating grout). When ``k_g = k_s``, ``\sigma = 0`` and the
  grout and ground are thermally indistinguishable.

The remaining groups ``\theta_i`` are geometric ratios that differ between the single- and
double-U configurations and are defined below.

## Single U-tube

With the geometric ratios

```math
\theta_1 = \frac{s}{2 r_b}, \qquad
\theta_2 = \frac{r_b}{r_o}, \qquad
\theta_3 = \frac{r_o}{s},
```

where ``\theta_1 = D/r_b`` is the dimensionless half-spacing (``D = s/2`` is the pipe offset
from the borehole centre).

### Zeroth order (line source)

Eq. 12 of Javed & Spitler (2017):

```math
R_b = \frac{1}{4\pi k_g}\left[\beta + \ln\!\frac{\theta_2}{2\,\theta_1\,(1-\theta_1^4)^{\sigma}}\right].
```

### First order

Eq. 13 of Javed & Spitler (2017) adds the first multipole correction. With ``b_1 = (1+\beta)/(1-\beta)``:

```math
R_b = \frac{1}{4\pi k_g}\left[
\beta + \ln\!\frac{\theta_2}{2\,\theta_1\,(1-\theta_1^4)^{\sigma}}
- \frac{\theta_3^2\left(1 - \dfrac{4\sigma\theta_1^4}{1-\theta_1^4}\right)^2}
       {b_1 + \theta_3^2\left(1 + \dfrac{16\sigma\theta_1^4}{(1-\theta_1^4)^2}\right)}
\right].
```

The first-order term is the recommended default (`order = 1`); the zeroth-order form is exposed
mainly for comparison with the classical line-source estimate.

## Double U-tube

For `nLoop = 2` (four pipes, the two loops sharing the borehole) the package uses the explicit
formulas of Claesson & Javed (2019). The zeroth-order borehole resistance is

```math
R_b = \frac{R_p^{\text{tot}}}{4}
+ \frac{1}{4\pi k_g}\left[
\ln\!\frac{r_b^4}{4\,r_o\,(s/2)^3}
+ \sigma \ln\!\frac{r_b^8}{r_b^8 - (s/2)^8}
\right],
```

with the ``R_p^{\text{tot}}/4`` term reflecting the four parallel pipes. The first-order form
adds a multipole correction built from

```math
\theta_1 = \frac{r_o^2}{4(s/2)^2}, \qquad
\theta_2 = \frac{(s/2)^2}{\bigl(r_b^8-(s/2)^8\bigr)^{1/4}}, \qquad
\theta_3 = \frac{r_b^2}{\bigl(r_b^8-(s/2)^8\bigr)^{1/4}},
```

(see the source for the full expression).

## Total internal resistance ``R_a``

``R_a`` is the resistance to heat exchange *between* the down-flowing and up-flowing legs — the
quantity that controls the thermal short-circuit. For the single U-tube, with
``\theta_1 = s/(2 r_b)`` and ``\theta_3 = r_o/s``, the zeroth-order form (Eq. 26 of Javed &
Spitler, 2017) is

```math
R_a = \frac{1}{\pi k_g}\left[\beta + \ln\!\frac{(1+\theta_1^2)^{\sigma}}{\theta_3\,(1-\theta_1^2)^{\sigma}}\right],
```

with a first-order correction analogous to ``R_b``. See
[`resistance_ULoop_total_internal`](@ref).

### Double-U pipe networks

For the double U-tube the two loops can be connected in two ways, selected with the `network`
keyword:

- `"diagonal"` (default) — the paired legs sit on the diagonal of the four-pipe arrangement
  (Eqs. 18–19 of Claesson & Javed, 2019);
- `"adjacent"` — the paired legs are neighbours (Eqs. 22–23).

!!! warning "Known limitation"
    The `nLoop = 2, order = 1, network = "adjacent"` branch can return an unphysical negative
    ``R_a`` for some geometries (flagged as a TODO in the source). Prefer the `"diagonal"`
    network, or fall back to `order = 0`, for adjacent-pair double U-tubes until this is
    resolved.

## Coaxial (concentric-tube) exchanger

Coaxial boreholes do not use the multipole network. Following Lamarche (2021), the cross-section
reduces to two resistances (implemented in
[`resistance_coaxial`](@ref)):

```math
R_{12} = \frac{1}{h_{in}\,\pi d_{ii}}
       + \frac{\ln(d_{io}/d_{ii})}{2\pi k_{p,in}}
       + \frac{1}{h_{ann}\,\pi d_{io}},
\qquad
R_1 = \frac{1}{h_{ann}\,\pi d_{oi}}
    + \frac{\ln(d_{oo}/d_{oi})}{2\pi k_{p,out}}
    + \frac{\ln(d_b/d_{oo})}{2\pi k_g},
```

(Eqs. 1–2 of Lamarche, 2021), where ``R_{12}`` links the center pipe to the annulus and ``R_1``
links the annulus fluid to the borehole wall. The convection coefficient in the center pipe uses
[`Nusselt`](@ref); the annulus uses [`Nusselt_annulus`](@ref). By Eq. 8, the (steady) borehole
resistance of a coaxial exchanger is simply ``R_b = R_1``.

### Effective resistance ``R_b^*``

With the groups (Eq. 7, ``\dot m c_f = V\rho_f c_f``)

```math
\gamma = \frac{H}{2\,\dot m c_f\,R_1}, \quad
R_a = \frac{4 R_1 R_{12}}{4 R_1 + R_{12}}, \quad
\xi = \sqrt{\frac{R_a}{4 R_1}}, \quad
\eta = \frac{\gamma}{\xi},
```

the two closed forms exposed by [`resistance_coaxial_effective`](@ref) are

```math
R_b^* = R_1\,\eta\coth\eta \quad (\text{UBW, Eq. 14}),
\qquad
R_b^* = R_1\left(1 + \frac{R_a}{R_{12}}\frac{\eta^2}{3}\right) \quad (\text{UHF, Eq. 31}).
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

## Recovering the grout-only resistance

Since ``R_b`` includes the fluid and pipe contributions, the grout-only resistance is the
remainder:

```math
R_g = R_b - R_p - R_f.
```

This is a useful sanity check — ``R_g`` must be positive.

## Functions on this page

```@docs
resistance_ULoop_borehole
resistance_ULoop_total_internal
resistance_coaxial
resistance_coaxial_effective
```
