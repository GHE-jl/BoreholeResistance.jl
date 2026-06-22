# Validation of thermal resistance functions for a single U-loop ground heat exchanger.
# Covers: Reynolds, Prandtl (both overloads), Nusselt (both overloads), resistance_fluid
# (both overloads), resistance_pipe, resistance_borehole_multipole (both overloads, order 0/1),
# resistance_total_internal_multipole (both overloads, order 0/1), and
# resistance_borehole_effective (all three overloads). Flow-rate sweep included.
# Run from package root: julia --project=script/ script/script_single_Uloop.jl

using BoreholeResistance

# --- Parameters (typical single U-loop borehole) ---
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

V_nom  = 15.0 / 1000 / 60         # 15 L/min in m³/s
V̇_nom = V_nom / (π * ri^2)       # Mean fluid speed [m/s]

# --- Point check at V = 15 L/min ---
println("=== Single U-loop — point check at V = 15 L/min ===")
Re_nom = Reynolds(V̇_nom, ri, ρf, μf)
Pr_nom = Prandtl(kf, cf, μf)

# Prandtl: mass-specific and volumetric overloads must agree
@assert Prandtl(kf, cf, μf) ≈ Prandtl(kf, cf * ρf, ρf, μf)  "Prandtl overloads must agree"

Nu_nom = Nusselt(Re_nom, Pr_nom, ri, ϵ)

# Nusselt: (Re, Pr, r, ϵ) and (V̇, r, ...) overloads must agree
@assert Nu_nom ≈ Nusselt(V̇_nom, ri, kf, cf, ρf, μf, ϵ)  "Nusselt overloads must agree"

Rp     = resistance_pipe(ro, ri, kp)
Rf_nom = resistance_fluid(V̇_nom, ri, kf, cf, ρf, μf, ϵ)

# resistance_fluid: (V̇, ...) and (Nu, r, kf) overloads must agree
@assert Rf_nom ≈ resistance_fluid(Nu_nom, ri, kf)  "resistance_fluid overloads must agree"

println("  Re = $(round(Int, Re_nom)),  Pr = $(round(Pr_nom, digits=3)),  Nu = $(round(Nu_nom, digits=3))")
println("  Rf = $(round(Rf_nom, digits=4)) mK/W")
println("  Rp = $(round(Rp, digits=4)) mK/W  (constant)")
println()

# --- Rb, Ra, Rbe for order 0 and 1 ---
println("=== Rb, Ra, Rbe at V = 15 L/min ===")
for order in [0, 1]
    Rb  = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
    Ra  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
    Rbe = resistance_borehole_effective(V_nom, H, cf, ρf, Rb, Ra)
    Rg  = Rb - Rp - Rf_nom

    # Long-form overloads (auto-compute Rf/Rp from pipe geometry) must agree
    @assert Rb ≈ resistance_borehole_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ; order=order)  "resistance_borehole_multipole overloads must agree"
    @assert Ra ≈ resistance_total_internal_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ; order=order)  "resistance_total_internal_multipole overloads must agree"

    @assert Rp     > 0
    @assert Rf_nom > 0
    @assert Rg     > 0  "Grout resistance must be positive"
    @assert Rb     > 0
    @assert Rbe   >= Rb  "Rbe must be ≥ Rb"

    println("  order $order:  Rb = $(round(Rb, digits=4))  Ra = $(round(Ra, digits=4))  Rbe = $(round(Rbe, digits=4))  Rg = $(round(Rg, digits=4)) mK/W")
end
println()

# resistance_borehole_effective: all three overloads must agree (order 1)
Rb1  = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
Ra1  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
Rbe1 = resistance_borehole_effective(V_nom, H, cf, ρf, Rb1, Ra1)
Rbe2 = resistance_borehole_effective(V_nom, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_nom)
Rbe3 = resistance_borehole_effective(V_nom, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ)
@assert Rbe1 ≈ Rbe2  "resistance_borehole_effective overloads 1 and 2 must agree"
@assert Rbe1 ≈ Rbe3  "resistance_borehole_effective overloads 1 and 3 must agree"

# --- Flow-rate sweep ---
println("=== Flow-rate sweep (single U-loop, order 1) ===")
println(rpad("V [L/min]", 12), " ", rpad("Re", 8), " ", rpad("Rf", 8), " ", rpad("Rp", 8), " ",
        rpad("Rg", 8), " ", rpad("Rb", 8), " ", "Rbe")
for Q_Lmin in [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]
    V   = Q_Lmin / 1000 / 60
    V̇  = V / (π * ri^2)
    Re  = Reynolds(V̇, ri, ρf, μf)
    Rf  = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ)
    Rb  = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf; order=1)
    Ra  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf; order=1)
    Rbe = resistance_borehole_effective(V, H, cf, ρf, Rb, Ra)
    Rg  = Rb - Rp - Rf
    println(rpad(string(round(Q_Lmin, digits=1)), 12), " ", rpad(string(round(Int, Re)), 8), " ",
            rpad(string(round(Rf, digits=4)), 8), " ", rpad(string(round(Rp, digits=4)), 8), " ",
            rpad(string(round(Rg, digits=4)), 8), " ", rpad(string(round(Rb, digits=4)), 8), " ",
            round(Rbe, digits=4))
end

println("\nAll assertions passed.")
