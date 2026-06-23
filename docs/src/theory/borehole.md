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
[`resistance_total_internal_multipole`](@ref).

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

## Recovering the grout-only resistance

Since ``R_b`` includes the fluid and pipe contributions, the grout-only resistance is the
remainder:

```math
R_g = R_b - R_p - R_f.
```

This is a useful sanity check — ``R_g`` must be positive.

## Functions on this page

```@docs
resistance_borehole_multipole
resistance_total_internal_multipole
```
