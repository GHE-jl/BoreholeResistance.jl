# References

The models implemented in this package are drawn from the following sources.

## Multipole method and borehole resistance

- **Hellström, G.** (1991). *Ground Heat Storage: Thermal Analyses of Duct Storage Systems.*
  Doctoral thesis, Lund University.
  The original derivation of the multipole method. The single-U borehole resistance
  corresponds to Eq. 8.36.

- **Javed, S., & Spitler, J.** (2017). Accuracy of borehole thermal resistance calculation
  methods for grouted single U-tube ground heat exchangers. *Applied Energy*, 187, 790–806.
  <https://doi.org/10.1016/j.apenergy.2016.11.079>
  Source of the explicit zeroth- and first-order formulas for ``R_b`` (Eqs. 12–13) and for the
  total internal resistance ``R_a`` (Eqs. 25–26) used for the single U-tube.

- **Claesson, J., & Javed, S.** (2019). Explicit multipole formulas and thermal network models
  for calculating thermal resistances of double U-pipe borehole heat exchangers.
  *Science and Technology for the Built Environment*, 25(8), 980–992.
  <https://doi.org/10.1080/23744731.2019.1620565>
  Source of the double U-tube formulas, including the diagonal (Eqs. 18–19) and adjacent
  (Eqs. 22–23) pipe-network variants.

- **Claesson, J., & Hellström, G.** (2011). Multipole method to calculate borehole thermal
  resistances in a borehole heat exchanger. *HVAC&R Research*, 17(6), 895–911.

## Effective resistance formulas

- **Javed, S., & Spitler, J. D.** (2016). Calculation of borehole thermal resistance.
  In S. J. Rees (Ed.), *Advances in Ground-Source Heat Pump Systems* (pp. 63–95).
  Woodhead Publishing. <https://doi.org/10.1016/B978-0-08-100311-4.00003-0>
  Source of the effective borehole resistance ``R_b^*`` under the UBW (Eqs. 3.68–3.70) and UHF
  (Eq. 3.67) boundary conditions.

## Convective heat transfer and fluid mechanics

- **Lamarche, L.** (2023). *Fundamentals of Geothermal Heat Pump Systems: Design and
  Application.* Springer Nature Switzerland.
  Source of the Nusselt-number treatment (Eqs. 2.42, 2.43b, 2.48, 2.49) and the convective
  resistance (Eqs. 2.32, 5.6).

- **Lamarche, L.** (2021). Analytic models and effective resistances for coaxial ground heat
  exchangers. *Geothermics*, 97, 102224.
  <https://doi.org/10.1016/j.geothermics.2021.102224>
  Source of the annulus Nusselt correlation (Eqs. 64a–64b).

- **Bergman, T. L., & Incropera, F. P.** (2011). *Fundamentals of Heat and Mass Transfer*
  (7th ed.). Wiley, New York.
  Cylindrical-shell conduction resistance used for ``R_p``.

## Thermophysical property data

- **The Engineering ToolBox.** Temperature-dependent data for the thermal conductivity, specific
  heat, density and dynamic viscosity of liquid water, to which the polynomial fits in
  [Water properties](@ref) were calibrated. <https://www.engineeringtoolbox.com>
