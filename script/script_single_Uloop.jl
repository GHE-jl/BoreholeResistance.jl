# Validation of the single U-loop thermal-resistance functions against Javed and Spitler (2017),
# "Accuracy of borehole thermal resistance calculation methods for grouted single U-tube ground
# heat exchangers", Applied Energy 187, 790-806.
#
# The zeroth- and first-order multipole borehole resistance (Eqs. 12-13) and the total internal
# resistance (Eqs. 25-26) are reproduced. The grout resistance Rg = Rb - Rp/N (N = 2, Eq. 3) is
# compared with Table 6 of the paper, which lists Rg for a fixed pipe resistance Rp.

using BoreholeResistance

# Table 6 geometry of Javed and Spitler (2017): pipe OD = 32 mm, ID = 28 mm; θ2 = rb/rpo = 3 gives
# borehole diameter 96 mm.
# Shank spacing is set through θ1 = s/(2·rb): close, moderate and wide (Configurations A, B, C).
# The pipe resistance Rp is prescribed directly (0.050 turbulent, 0.185 laminar) — it is passed as
# the code's `Rp` while the fluid resistance is set to 0 so that Rp + Rf equals the paper's Rp.
rpo = 0.032 / 2         # Pipe outer radius [m]
rb  = 3 * rpo           # θ2 = 3 → borehole radius [m]
ks  = 3.0               # Ground thermal conductivity [W/mK]

# θ1 values for θ2 = 3 (Table 2 of the paper): close / moderate / wide
spacings = [("Close", 0.333), ("Moderate", 0.445), ("Wide", 0.667)]

# Reference grout resistances Rg [mK/W] from Table 6, keyed by (kg, Rp, spacing). Results are the
# zeroth- and first-order multipole values.
ref_Rg = Dict(
    (0.6, 0.050, "Close")    => [0.19839, 0.17732],
    (0.6, 0.050, "Moderate") => [0.15783, 0.14297],
    (0.6, 0.050, "Wide")     => [0.08810, 0.07090],
    (0.6, 0.185, "Close")    => [0.19839, 0.19231],
    (0.6, 0.185, "Moderate") => [0.15783, 0.15380],
    (0.6, 0.185, "Wide")     => [0.08810, 0.08394],
    (2.4, 0.050, "Close")    => [0.04983, 0.04869],
    (2.4, 0.050, "Moderate") => [0.04019, 0.03953],
    (2.4, 0.050, "Wide")     => [0.02608, 0.02572],
    (2.4, 0.185, "Close")    => [0.04983, 0.05430],
    (2.4, 0.185, "Moderate") => [0.04019, 0.04262],
    (2.4, 0.185, "Wide")     => [0.02608, 0.02730],
)

println("=== Single U-loop — grout resistance vs Javed and Spitler (2017), Table 6 ===")
println(rpad("kg", 5), " ", rpad("Rp", 6), " ", rpad("spacing", 9), " ",
        rpad("Rg(0)", 8), " ", rpad("Rg(0) paper", 9), " ", rpad("err0%", 7), " ",
        rpad("Rg(1)", 8), " ", rpad("Rg(1) paper", 9), " ", rpad("err1%", 7))
for kg in [0.6, 2.4], Rp in [0.050, 0.185], (name, θ1) in spacings
    s = 2 * rb * θ1                       # shank spacing from θ1 = s/(2·rb)
    Rg0 = resistance_ULoop_borehole(s, rb, rpo, ks, kg, Rp, 0.0; order=0) - Rp / 2
    Rg1 = resistance_ULoop_borehole(s, rb, rpo, ks, kg, Rp, 0.0; order=1) - Rp / 2
    Rg_ref = ref_Rg[(kg, Rp, name)]
    err0 = 100 * abs(Rg0 - Rg_ref[1]) / Rg_ref[1]
    err1 = 100 * abs(Rg1 - Rg_ref[2]) / Rg_ref[2]
    println(rpad(kg, 5), " ", rpad(Rp, 6), " ", rpad(name, 9), " ",
            rpad(round(Rg0, digits=5), 8), " ", rpad(round(Rg_ref[1], digits=5), 9), " ", rpad(round(err0, digits=2), 7), " ",
            rpad(round(Rg1, digits=5), 8), " ", rpad(round(Rg_ref[2], digits=5), 9), " ", rpad(round(err1, digits=2), 7))
    # The implemented first-order formula reproduces the paper's first-order column (residual is
    # only from the 3-digit rounding of θ1 in the paper)
    @assert err1 < 1.0  "Rg(1) must match Table 6 first-order column ($name, kg=$kg, Rp=$Rp)"
end
println()

# Realistic HDPE single U-loop with water properties
H  = 100.0          # Borehole length [m]
rb = 0.075          # Borehole radius [m]
ro = 0.020          # Pipe outer radius [m]
ri = 0.0164         # Pipe inner radius [m]
s  = 0.08           # Shank spacing (centre-to-centre) [m]
T0 = 10.0           # Reference fluid temperature [°C]
ks = 2.0            # Ground thermal conductivity [W/mK]
kg = 1.0            # Grout thermal conductivity [W/mK]
kp = 0.4            # Pipe thermal conductivity [W/mK]
ϵ  = 5e-6           # Pipe roughness [m]

kf = water_k(T0)
cf = water_cp(T0)
ρf = water_ρ(T0)
μf = water_μ(T0)

V_nom  = 15.0 / 1000 / 60           # 15 L/min in m³/s
V̇_nom = V_nom / (π * ri^2)          # Mean fluid speed [m/s]

Re_nom = Reynolds(V̇_nom, ri, ρf, μf)
Pr_nom = Prandtl(kf, cf, μf)
Nu_nom = Nusselt(Re_nom, Pr_nom, ri, ϵ)

# Overload agreement checks
@assert Prandtl(kf, cf, μf) ≈ Prandtl(kf, cf * ρf, ρf, μf)  "Prandtl overloads must agree"
@assert Nu_nom ≈ Nusselt(V̇_nom, ri, kf, cf, ρf, μf, ϵ)      "Nusselt overloads must agree"

Rp = resistance_pipe(ro, ri, kp)
Rf_nom = resistance_fluid(V̇_nom, ri, kf, cf, ρf, μf, ϵ)
@assert Rf_nom ≈ resistance_fluid(Nu_nom, ri, kf)           "resistance_fluid overloads must agree"

println("=== Single U-loop — HDPE example at V = 15 L/min ===")
println("  Re = $(round(Int, Re_nom)),  Pr = $(round(Pr_nom, digits=3)),  Nu = $(round(Nu_nom, digits=3))")
println("  Rf = $(round(Rf_nom, digits=4)) mK/W,  Rp = $(round(Rp, digits=4)) mK/W")
println()

println("=== Rb, Ra, Rbe at V = 15 L/min ===")
for order in [0, 1]
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
    Ra  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
    Rbe = resistance_ULoop_effective(V_nom, H, cf, ρf, Rb, Ra)
    Rg  = Rb - (Rp + Rf_nom) / 2            # grout resistance, N = 2 (Eq. 3)

    @assert Rb ≈ resistance_ULoop_borehole(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ;
        order=order)  "resistance_ULoop_borehole overloads must agree"
    @assert Ra ≈ resistance_ULoop_total_internal(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf,
        μf, ϵ; order=order)  "resistance_ULoop_total_internal overloads must agree"
    @assert Rg  > 0
    @assert Rbe >= Rb  "Rbe must be ≥ Rb"

    println("  order $order:  Rb = $(round(Rb, digits=4))  Ra = $(round(Ra, digits=4))  " *
            "Rbe = $(round(Rbe, digits=4))  Rg = $(round(Rg, digits=4)) mK/W")
end

# resistance_ULoop_effective — all three overloads must agree (order 1, mean model)
Rb1  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
Ra1  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
Rbe1 = resistance_ULoop_effective(V_nom, H, cf, ρf, Rb1, Ra1)
Rbe2 = resistance_ULoop_effective(V_nom, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_nom)
Rbe3 = resistance_ULoop_effective(V_nom, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ)
@assert Rbe1 ≈ Rbe2  "resistance_ULoop_effective overloads 1 and 2 must agree"
@assert Rbe1 ≈ Rbe3  "resistance_ULoop_effective overloads 1 and 3 must agree"
println()

# Flow-rate sweep
println("=== Flow-rate sweep (single U-loop, order 1) ===")
println(rpad("V [L/min]", 12), " ", rpad("Re", 8), " ", rpad("Rf", 8), " ", rpad("Rp", 8), " ",
        rpad("Rg", 8), " ", rpad("Rb", 8), " ", "Rbe")
for Q_Lmin in [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]
    V   = Q_Lmin / 1000 / 60
    V̇  = V / (π * ri^2)
    Re  = Reynolds(V̇, ri, ρf, μf)
    Rf  = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ)
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; order=1)
    Ra  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf; order=1)
    Rbe = resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra)
    Rg  = Rb - (Rp + Rf) / 2
    println(rpad(string(round(Q_Lmin, digits=1)), 12), " ", rpad(string(round(Int, Re)), 8), " ",
            rpad(string(round(Rf, digits=4)), 8), " ", rpad(string(round(Rp, digits=4)), 8), " ",
            rpad(string(round(Rg, digits=4)), 8), " ", rpad(string(round(Rb, digits=4)), 8), " ",
            round(Rbe, digits=4))
end
