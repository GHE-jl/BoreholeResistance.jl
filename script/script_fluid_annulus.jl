# Validation of dimensionless number and friction factor functions for annulus flow
# (coaxial borehole geometry). Covers: Reynolds, Prandtl (both overloads), Nusselt
# (both overloads), Nusselt_annulus (both overloads), and friction factor comparison
# between the iterative Colebrook-White and the explicit Tkachenko-Mileikovskyi correlation.

using BoreholeResistance

# Parameters (coaxial geometry: fluid in annulus between borehole and inner pipe)
rb = 0.075       # Borehole radius = outer annulus radius [m]
ro = 0.030       # Inner pipe outer radius = inner annulus radius [m]
T0 = 10.0        # Reference fluid temperature [°C]
ϵ  = 5e-6        # Pipe roughness [m]

kf, cf, ρf, μf = fluid_property(T0, :water)

r_h = rb - ro                   # Annulus hydraulic radius = D_h / 2 [m]
A_annulus = π * (rb^2 - ro^2)   # Annulus cross-sectional area [m²]

# Point check at V = 15 L/min
V_nom  = 15.0 / 1000 / 60
V̇_nom = V_nom / A_annulus

Re_nom = Reynolds(V̇_nom, r_h, ρf, μf)
Pr_nom = Prandtl(kf, cf, μf)

# Prandtl: mass-specific and volumetric overloads must agree
@assert Prandtl(kf, cf, μf) ≈ Prandtl(kf, cf * ρf, ρf, μf)  "Prandtl overloads must agree"

println("=== Point check at V = 15 L/min ===")
println("  Re = $(round(Int, Re_nom)),  Pr = $(round(Pr_nom, digits=3))")
println()

# Friction factors: Colebrook-White vs Tkachenko-Mileikovskyi
println("=== Friction factor comparison (Colebrook-White vs Tkachenko-Mileikovskyi) ===")
println(rpad("Re", 10), " ", rpad("f_CW", 12), " ", rpad("f_TM", 12), " ", "rel. diff")
for Re_ff in [1000.0, 2300.0, 4000.0, 1e4, 5e4, 1e5, 5e5]
    f_CW = friction_factor_Colebrook_White(Re_ff, r_h, ϵ)
    f_TM = friction_factor_Tkachenko_Mileikovskyi(Re_ff, r_h, ϵ)
    rel_diff = f_CW > 0 ? abs(f_CW - f_TM) / f_CW : 0.0
    println(rpad(string(round(Int, Re_ff)), 10), " ",
            rpad(string(round(f_CW, digits=6)), 12), " ",
            rpad(string(round(f_TM, digits=6)), 12), " ",
            round(rel_diff, digits=5))
    if Re_ff >= 4000
        @assert rel_diff < 0.01  "TM and CW must agree within 1% for turbulent flow (Re=$Re_ff)"
    end
end
println()

# Nusselt: pipe formula vs annulus-specific correlation
Nu_pipe = Nusselt(Re_nom, Pr_nom, r_h, ϵ)
Nu_ann  = Nusselt_annulus(Re_nom, Pr_nom, rb, ro, ϵ, ϵ)

# Nusselt overloads must agree
@assert Nu_pipe ≈ Nusselt(V̇_nom, r_h, kf, cf, ρf, μf, ϵ)  "Nusselt overloads must agree"
@assert Nu_ann  ≈ Nusselt_annulus(V̇_nom, rb, ro, kf, cf, ρf, μf, ϵ, ϵ)  "Nusselt_annulus overloads must agree"

println("=== Nusselt at V = 15 L/min ===")
println("  Nusselt (pipe formula, D_h)  = $(round(Nu_pipe, digits=4))")
println("  Nusselt_annulus (Re, Pr)     = $(round(Nu_ann, digits=4))")
println()

# Flow-rate sweep
println("=== Flow-rate sweep ===")
println(rpad("V [L/min]", 12), " ", rpad("Re", 8), " ", rpad("Nu_pipe", 14), " ",
        rpad("Nu_annulus", 14), " ", "ratio")
for Q_Lmin in [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]
    V  = Q_Lmin / 1000 / 60
    V̇ = V / A_annulus
    Re = Reynolds(V̇, r_h, ρf, μf)
    Pr = Prandtl(kf, cf, μf)
    Nu1 = Nusselt(Re, Pr, r_h, ϵ)
    Nu2 = Nusselt_annulus(Re, Pr, rb, ro, ϵ, ϵ)
    println(rpad(string(round(Q_Lmin, digits=1)), 12), " ", rpad(string(round(Int, Re)), 8), " ",
            rpad(string(round(Nu1, digits=4)), 14), " ", rpad(string(round(Nu2, digits=4)), 14), " ",
            round(Nu2/Nu1, digits=4))
end
