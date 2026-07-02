# Validation of thermal resistance functions for a double U-loop ground heat exchanger.
# Covers: resistance_ULoop_borehole (both overloads, nLoop=2, order 0/1),
# resistance_ULoop_total_internal (both overloads, diagonal and adjacent networks),
# resistance_ULoop_effective, and single-vs-double Rb comparison.
# Run from package root: julia --project=script/ script/script_double_Uloop.jl

using BoreholeResistance

# --- Parameters ---
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

V_nom  = 15.0 / 1000 / 60
V̇_nom = V_nom / (π * ri^2)

Rp = resistance_pipe(ro, ri, kp)
Rf = resistance_fluid(V̇_nom, ri, kf, cf, ρf, μf, ϵ)

# --- Point check at V = 15 L/min ---
println("=== Double U-loop — point check at V = 15 L/min ===")
println("  Rf = $(round(Rf, digits=4)) mK/W,  Rp = $(round(Rp, digits=4)) mK/W")
println()

for order in [0, 1]
    Rb      = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=2, order=order)
    Ra_diag = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf;
                  nLoop=2, order=order, network="diagonal")
    Ra_adj  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf;
                  nLoop=2, order=order, network="adjacent")

    # Long-form overloads must agree with short-form
    @assert Rb ≈ resistance_ULoop_borehole(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ;
        nLoop=2, order=order)  "resistance_ULoop_borehole overloads must agree (order $order)"
    @assert Ra_diag ≈ resistance_ULoop_total_internal(V_nom, s, rb, ro, ri, ks, kg, kp, kf, cf,
        ρf, μf, ϵ; nLoop=2, order=order, network="diagonal")  "resistance_ULoop_total_internal
        overloads must agree"

    @assert Rb > 0
    @assert Ra_diag > 0
    if order == 0; @assert Ra_adj > 0; end  # order-1 adjacent can be negative for typical params

    println("  order $order:  Rb = $(round(Rb, digits=4)) mK/W")
    println("            Ra (diagonal) = $(round(Ra_diag, digits=4)) mK/W")
    print("            Ra (adjacent) = $(round(Ra_adj, digits=4)) mK/W")
    if Ra_adj <= 0
        print("  ⚠ negative — expected for order 1 with these parameters")
    end
    println()
    println()
end

# --- Rbe for double U-loop (diagonal, order 1) ---
println("=== Rbe for double U-loop (order 1, diagonal) ===")
Rb_d  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=2, order=1)
Ra_d  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf;
            nLoop=2, order=1, network="diagonal")
Rbe_d = resistance_ULoop_effective(V_nom, H, cf, ρf, Rb_d, Ra_d)
@assert Rbe_d >= Rb_d  "Rbe must be ≥ Rb"
println("  Rb = $(round(Rb_d, digits=4)),  Ra = $(round(Ra_d, digits=4)),
    Rbe = $(round(Rbe_d, digits=4)) mK/W")
println()

# --- Single vs double comparison (order 1) ---
println("=== Single vs double U-loop comparison (order 1) ===")
Rb1 = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=1, order=1)
Rb2 = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=2, order=1)
println("  Rb single = $(round(Rb1, digits=4)) mK/W")
println("  Rb double = $(round(Rb2, digits=4)) mK/W")
@assert Rb1 > 0 && Rb2 > 0
println()

# --- Flow-rate sweep (diagonal, order 1) ---
println("=== Flow-rate sweep (double U-loop, diagonal, order 1) ===")
println(rpad("V [L/min]", 12), " ", rpad("Re", 8), " ", rpad("Rf", 8), " ", rpad("Rb", 8), " ",
    "Ra_diag")
for Q_Lmin in [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]
    local V, V̇, Re, Rf_i, Rb, Ra
    V    = Q_Lmin / 1000 / 60
    V̇   = V / (π * ri^2)
    Re   = Reynolds(V̇, ri, ρf, μf)
    Rf_i = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ)
    Rb   = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_i; nLoop=2, order=1)
    Ra   = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_i;
               nLoop=2, order=1, network="diagonal")
    println(rpad(string(round(Q_Lmin, digits=1)), 12), " ", rpad(string(round(Int, Re)), 8), " ",
            rpad(string(round(Rf_i, digits=4)), 8), " ", rpad(string(round(Rb, digits=4)), 8), " ",
            round(Ra, digits=4))
end

println("\nAll assertions passed.")
