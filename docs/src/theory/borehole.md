# Borehole (grout) resistance

The borehole resistance ``R_b`` combines the fluid, pipe and grout contributions into the
single resistance from the fluid to the borehole wall. Because the pipes sit off-centre and
interact, the grout step is solved with the **multipole method** of Hellström (1991), in the
explicit form given by Javed & Spitler (2017) for single U-tubes and Claesson & Javed (2019)
for double U-tubes. The same machinery also yields the **total internal resistance** ``R_a``
between the two legs. Both are implemented in
[`resistance_borehole.jl`](https://github.com/GHE-jl/BoreholeResistance.jl/blob/main/src/resistance_borehole.jl).

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
R_b^{(0)} = \frac{1}{4\pi k_g}\left[\beta + \ln\!\frac{\theta_2}{2\,\theta_1\,(1-\theta_1^4)^{\sigma}}\right].
```

### First order

Eq. 13 of Javed & Spitler (2017) adds the first multipole correction. With ``b_1 = (1+\beta)/(1-\beta)``:

```math
R_b^{(1)} = \frac{1}{4\pi k_g}\left[
\beta + \ln\!\frac{\theta_2}{2\,\theta_1\,(1-\theta_1^4)^{\sigma}}
- \frac{\theta_3^2\left(1 - \dfrac{4\sigma\theta_1^4}{1-\theta_1^4}\right)^2}
       {b_1 + \theta_3^2\left(1 + \dfrac{16\sigma\theta_1^4}{(1-\theta_1^4)^2}\right)}
\right].
```

The first-order term is the recommended default (`order = 1`); the zeroth-order form is exposed
mainly for comparison with the classical line-source estimate.

## Double U-tube

For `nLoop = 2` (four pipes, the two loops sharing the borehole) the package uses the explicit formulas of Claesson & Javed (2019).

### Zeroth order

 The zeroth-order borehole resistance (Eq. 13) is

```math
R_b^{(0)} = \frac{R_p^{\text{tot}}}{4}
+ \frac{1}{8\pi k_g}\left[
\ln\!\frac{r_b^4}{4\,r_o\,(s/2)^3}
+ \sigma \ln\!\frac{r_b^8}{r_b^8 - (s/2)^8}
\right],
```

with the ``R_p^{\text{tot}}/4`` term reflecting the four parallel pipes. The package
parameterises the four-pipe geometry through the shank spacing ``s``, the distance between diagonally opposite pipes.

### First order

The first-order form adds a multipole correction built from

```math
\theta_1 = \frac{r_o^2}{4(s/2)^2}, \qquad
\theta_2 = \frac{(s/2)^2}{\bigl(r_b^8-(s/2)^8\bigr)^{1/4}}, \qquad
\theta_3 = \frac{r_b^2}{\bigl(r_b^8-(s/2)^8\bigr)^{1/4}},
```

Then, the explicit formula adds to the line source model with

```math
R_b^{(1)} = R_b^{(0)} - \frac{1}{8\pi k_g}\cdot\frac{b_1 \theta_1 \left(3-8\sigma \theta_2^4\right)^2}{1+b_1 \theta_1 \left(5+64\sigma \theta_2^4 \theta_3^4 \right)}.
```

## Total internal resistance ``R_a``

``R_a`` is the resistance to heat exchange *between* the down-flowing and up-flowing legs, the quantity that controls the thermal short-circuit.

### Single U-pipe network

For the single U-tube, with ``\theta_1 = s/(2 r_b)`` and ``\theta_3 = r_o/s``, the zeroth-order form (Eq. 25 of Javed & Spitler, 2017) is

```math
R_a^{(0)} = \frac{1}{\pi k_g}\left[\beta + \ln\left(\frac{(1+\theta_1^2)^{\sigma}}{\theta_3\,(1-\theta_1^2)^{\sigma}}\right)\right],
```

with a first-order correction analogous to ``R_b`` (Eq. 26 of Javed & Spitler, 2017)

```math
R_a^{(1)} = \frac{1}{\pi k_g}\left[\beta + \ln\left(\frac{(1+\theta_1^2)^{\sigma}}{\theta_3\,(1-\theta_1^2)^{\sigma}}\right) - \frac{\theta_3^2\left(1-\theta_1^4+4\sigma\theta_1^2\right)^2}{b_1\left(1-\theta_1^4\right)^2 - \theta_3^2\left(1-\theta_1^4\right)^2+8\sigma \theta_1^2 \theta_3^2 \left(1+\theta_1^4\right)} \right].
```

### Double U-pipe networks

For the double U-tube the two loops can be connected in two ways, selected with the `network`
keyword:

- `"diagonal"` (default) — the paired legs sit on the diagonal of the four-pipe arrangement (Eqs. 18–19 of Claesson & Javed, 2019);
- `"adjacent"` — the paired legs are neighbours (Eqs. 22–23).

For the `"diagonal"` network, the zeroth order is (Eq. 18 of Claesson & Javed, 2019)

```math
R_{a,d}^{(0)} = 2R_p^{\text{tot}} + \frac{1}{\pi k_g}\left(\ln\left(\frac{(s/2)}{r_o}\right) + \sigma \ln\left(\frac{r_b^4+(s/2)^4}{r_b^4 - (s/2)^4}\right)\right),
```

while the first order is (Eq. 19 of Claesson & Javed, 2019)

```math
R_{a,d}^{(1)} = R_{a,d}^{(0)}-\frac{1}{\pi k_g}\cdot\frac{b_1\theta_1\left(1+8\sigma\theta_2^2\theta_3^2\right)^2}{1-b_1\theta_1\left(3-32\sigma\left(\theta_2^2\theta_3^6+\theta_2^6\theta_3^2\right)\right)}.
```

For the `"adjacent"` network, the zeroth order is (Eq. 22 of Claesson & Javed, 2019)

```math
R_{a,a}^{(0)} = 2R_p^{\text{tot}} + \frac{1}{\pi k_g}\left(\ln\left(\frac{s}{r_o}\right) + \sigma \ln\left(\frac{r_b^4+(s/2)^4}{r_b^4 - (s/2)^4}\right)\right),
```

while the first order is (Eq. 23 of Claesson & Javed, 2019)

```math
R_{a,a}^{(1)} = R_{a,a}^{(0)}+\frac{b_1\theta_1}{2\pi k_g}\cdot\frac{V_2^2M_{11}-2V_1V_2M_{21}-V_1^2M_{22}}{M_{11}M_{22}+M_{21}^2},
```

with the following variables

```math
V_1 = 1-8\sigma\theta_2^3\theta_3\\
V_2 = 3+8\sigma\theta_2\theta_3^3\\
M_{11} = 1+16b_1\sigma\theta_1\left(3\theta_2^3\theta_3^5+\theta_2^7\theta_3\right)\\
M_{12} = -M_{21}\\
M_{21} = b_1 \theta_1\\
M_{22} = -1-16b_1\sigma\theta_1\left(\theta_2\theta_3^7+3\theta_2^5\theta_3^3\right)
```

## Coaxial (concentric-tube) exchanger

Coaxial boreholes do not use the multipole network. Following Lamarche (2021), the cross-section
reduces to two resistances (implemented in
[`resistance_coaxial`](@ref)):

```math
R_{12} = \frac{1}{2\pi h_{in} r_{ii}}
       + \frac{\ln(r_{io}/r_{ii})}{2\pi k_{p,in}}
       + \frac{1}{2\pi h_{ann}r_{io}},
\qquad
R_1 = \frac{1}{2\pi h_{ann}r_{oi}}
    + \frac{\ln(r_{oo}/r_{oi})}{2\pi k_{p,out}}
    + \frac{\ln(r_b/r_{oo})}{2\pi k_g},
```

(Eqs. 1–2 of Lamarche, 2021), where ``R_{12}`` links the center pipe to the annulus fluid and ``R_1``
links the annulus fluid to the borehole wall. The convection coefficient in the center pipe uses
[`Nusselt`](@ref), and the annulus uses [`Nusselt_annulus`](@ref). By Eq. 8, the (steady) borehole
resistance of a coaxial exchanger is simply

```math
R_b = R_1,
```
and the total internal resistance is

```math
R_a = \frac{4 R_1 R_{12}}{4 R_1 + R_{12}}.
```

## Recovering the grout-only resistance

Since ``R_b`` includes the fluid and pipe contributions, the grout-only resistance is the
remainder. The ``N`` parallel pipes (``N = 2`` for a single U-tube, ``N = 4`` for a double)
contribute the combined fluid-and-pipe resistance ``R_p^{\text{tot}}/N``, so

```math
R_g = R_b - \frac{R_p + R_f}{N}
```

(Eq. 3 of Javed & Spitler, 2017). This is a useful sanity check — ``R_g`` must be positive.

## Functions on this page

```@docs
resistance_ULoop_borehole
resistance_ULoop_total_internal
resistance_coaxial
```
