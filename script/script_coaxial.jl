# Validation of thermal resistance functions for a coaxial (concentric-tube) ground heat
# exchanger. Covers: resistance_coaxial (both overloads: convection coefficients and V + fluid
# properties) and resistance_coaxial_effective (both overloads, UHF/UBW/mean/UHF_gradient models).
# The geometry is "Borehole 1" of Lamarche (2021, Table 1), so the results can be compared with
# the published R'12, R'1 (Section 2.3) and effective resistances (Table 2). The "UHF_gradient"
# section below uses "Borehole 3b" (Table 1/2) to compare against Table 3.
# Run from package root: julia --project=script/ script/script_coaxial.jl

using BoreholeResistance

# Parameters — Borehole 1 of Lamarche (2021), Table 1 (diameters → radii)
H   = 150.0          # Borehole length [m]
rb  = 0.1524 / 2     # Borehole radius [m]
roo = 0.1143 / 2     # Outer pipe outer radius [m]
roi = 0.0973 / 2     # Outer pipe inner radius [m]
rio = 0.0603 / 2     # Inner pipe outer radius [m]
rii = 0.0494 / 2     # Inner pipe inner radius [m]
kpi = 0.4            # Inner pipe thermal conductivity [W/mK]
kpo = 0.4            # Outer pipe thermal conductivity [W/mK]
kg  = 1.7            # Grout thermal conductivity [W/mK]
ks  = 2.5            # Ground thermal conductivity [W/mK] (not needed for Rb, kept for reference)
ϵ   = 5e-6           # Pipe roughness [m]
T0  = 10.0           # Reference fluid temperature [°C]
ṁ  = 3.25           # Mass flow rate [kg/s]

kf = water_k(T0)
cf = water_cp(T0)
ρf = water_ρ(T0)
μf = water_μ(T0)

V_nom = ṁ / ρf       # Volumetric flow rate [m³/s]

# --- Point check against the published resistances (Lamarche 2021, h given in Section 2.3) ---
# The paper reports hin = 5170, hann = 2150 W/m²K, giving R'12 = 0.084 and R'1 = 0.092 mK/W.
hin_paper  = 5170.0
hann_paper = 2150.0
R1_paper, R12_paper = resistance_coaxial(rii, rio, roi, roo, rb, kg, kpi, kpo, hin_paper,
    hann_paper)

println("=== Coaxial resistances — Borehole 1 (Lamarche 2021, published h) ===")
println("  R12 = $(round(R12_paper, digits=4)) mK/W   (paper: 0.084)")
println("  R1  = $(round(R1_paper, digits=4)) mK/W   (paper: 0.092)")
@assert isapprox(R12_paper, 0.084; atol=0.002)  "R12 must match the published value"
@assert isapprox(R1_paper, 0.092; atol=0.002)   "R1 must match the published value"
println()

# --- Resistances computed from the built-in Nusselt correlations ---
R1, R12 = resistance_coaxial(V_nom, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ)

println("=== Coaxial resistances — built-in Nusselt correlations at ṁ = $ṁ kg/s ===")
println("  R12 = $(round(R12, digits=4)) mK/W")
println("  R1  = $(round(R1, digits=4)) mK/W")
@assert R12 > 0
@assert R1  > 0
println()

# --- Effective borehole resistance (UHF / UBW / mean) ---
# Section 4 (Comsol comparison, Table 2) neglects the convection resistances, giving
# R'1 = 0.091 and R'12 = 0.0798 mK/W. Reproduce it with h → ∞ to compare against Table 2.
R1_t2, R12_t2 = resistance_coaxial(rii, rio, roi, roo, rb, kg, kpi, kpo, 1e12, 1e12)
println("=== Effective borehole resistance — Borehole 1 (Table 2, convection neglected) ===")
println("  R1 = $(round(R1_t2, digits=4)) (paper 0.091), R12 = $(round(R12_t2, digits=4)) (paper 0.0798)")
for model in ["UHF", "UBW", "mean"]
    Rbe = resistance_coaxial_effective(V_nom, H, cf, ρf, R1_t2, R12_t2; model=model)
    println("  model $(rpad(model, 5)): Rbe = $(round(Rbe, digits=4)) mK/W")
    @assert Rbe >= R1_t2  "Short-circuiting can only increase the effective resistance"
end
# Paper (Table 2, Borehole 1): Rb* ≈ 0.0915 (Eq. 31, UHF), 0.0916 (Eq. 14, UBW)
Rbe_UHF = resistance_coaxial_effective(V_nom, H, cf, ρf, R1_t2, R12_t2; model="UHF")
Rbe_UBW = resistance_coaxial_effective(V_nom, H, cf, ρf, R1_t2, R12_t2; model="UBW")
@assert isapprox(Rbe_UHF, 0.0915; atol=0.001)  "UHF effective resistance must match Table 2"
@assert isapprox(Rbe_UBW, 0.0916; atol=0.001)  "UBW effective resistance must match Table 2"

# Both overloads must agree (full form recomputes R1/R12 from geometry)
Rbe_full = resistance_coaxial_effective(V_nom, H, rii, rio, roi, roo, rb, kg, kpi, kpo,
    kf, cf, ρf, μf, ϵ; model="UHF")
Rbe_short = resistance_coaxial_effective(V_nom, H, cf, ρf, R1, R12; model="UHF")
@assert Rbe_full ≈ Rbe_short  "resistance_coaxial_effective overloads must agree"
println()

# --- Flow-rate sweep (built-in correlations, UHF model) ---
A_ann = π * (roi^2 - rio^2)       # Annulus cross-sectional area [m²]
println("=== Flow-rate sweep (coaxial, UHF model) ===")
println(rpad("ṁ [kg/s]", 12), " ", rpad("Re_ann", 10), " ", rpad("Re_in", 10), " ",
        rpad("R12", 8), " ", rpad("R1", 8), " ", "Rbe")
for ṁ_i in [0.4, 0.8, 1.6, 3.25, 5.0]
    V   = ṁ_i / ρf
    Re_ann = Reynolds(V / A_ann, roi - rio, ρf, μf)
    Re_in  = Reynolds(V / (π * rii^2), rii, ρf, μf)
    R1_i, R12_i = resistance_coaxial(V, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ)
    Rbe = resistance_coaxial_effective(V, H, cf, ρf, R1_i, R12_i; model="UHF")
    println(rpad(string(round(ṁ_i, digits=2)), 12), " ",
            rpad(string(round(Int, Re_ann)), 10), " ", rpad(string(round(Int, Re_in)), 10), " ",
            rpad(string(round(R12_i, digits=4)), 8), " ", rpad(string(round(R1_i, digits=4)), 8),
            " ", round(Rbe, digits=4))
end

# --- UHF_gradient model — Borehole 3b, Table 3 (linearly-varying far-field temperature) ---
# Table 3 case 1 ("Heat injection, Annulus-in") reports Rb* = 0.0420 mK/W via Eq. 61, versus
# 0.0356 mK/W for the plain UHF model (Eq. 31, gradient-independent by construction). This model
# needs no extra numeric input, but only applies to the flow-direction/heat-mode pairing that is
# "unfavorable" for the far-field gradient's sign — see the resistance_coaxial_effective docstring.
H3b  = 200.0
R1_3b, R12_3b = 0.0209, 0.0798      # Table 2 (Borehole 3b) and Section 4 (Borehole 3, R'12)
ṁ_3b = 0.8
V_3b = ṁ_3b / ρf

Rbe_UHF_3b  = resistance_coaxial_effective(V_3b, H3b, cf, ρf, R1_3b, R12_3b; model="UHF")
Rbe_grad_3b = resistance_coaxial_effective(V_3b, H3b, cf, ρf, R1_3b, R12_3b; model="UHF_gradient")
println("=== UHF vs UHF_gradient — Borehole 3b (Table 3, case 1: heat injection, annulus-in) ===")
println("  R1 = $R1_3b, R12 = $R12_3b mK/W")
println("  model UHF         : Rbe = $(round(Rbe_UHF_3b, digits=4)) mK/W   (paper Eq. 31: 0.0356)")
println("  model UHF_gradient: Rbe = $(round(Rbe_grad_3b, digits=4)) mK/W   (paper Eq. 61: 0.0420)")
@assert isapprox(Rbe_UHF_3b, 0.0356; atol=0.001)   "UHF effective resistance must match Table 3"
@assert isapprox(Rbe_grad_3b, 0.0420; atol=0.001)  "UHF_gradient effective resistance must match Table 3"
@assert Rbe_grad_3b > Rbe_UHF_3b  "the gradient correction amplifies the short-circuit resistance"
println()

println("\nAll assertions passed.")
