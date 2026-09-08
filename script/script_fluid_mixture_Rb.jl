# Synthetic example: effect of the heat-carrier fluid on the first-order multipole borehole
# resistance Rb of a single U-loop. Geometry and flow rate are held fixed; only the fluid (type and
# antifreeze concentration) changes, at a temperature (0°C) where pure water is no longer usable and
# an antifreeze mixture is required.

using BoreholeResistance

# Single U-loop HDPE geometry (same as script_single_Uloop.jl)
rb = 0.075          # Borehole radius [m]
ro = 0.020          # Pipe outer radius [m]
ri = 0.0164         # Pipe inner radius [m]
s  = 0.08           # Shank spacing (centre-to-centre) [m]
ks = 2.0            # Ground thermal conductivity [W/mK]
kg = 1.0            # Grout thermal conductivity [W/mK]
kp = 0.4            # Pipe thermal conductivity [W/mK]
ϵ  = 5e-6           # Pipe roughness [m]

T0 = 0.0                            # Reference fluid temperature [°C]
V  = 15.0 / 1000 / 60               # Flow rate, 15 L/min in m³/s
V̇  = V / (π * ri^2)                 # Mean fluid speed [m/s]

Rp = resistance_pipe(ro, ri, kp)

function Rb_for(fluid, percentage)
    kf, cf, ρf, μf = fluid == :water ? fluid_property(T0, :water) :
                                        fluid_property(T0, fluid; percentage=percentage)
    Re = Reynolds(V̇, ri, ρf, μf)
    Rf = resistance_fluid(V̇, ri, kf, cf, ρf, μf, ϵ)
    Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; order=1)
    return Re, Rf, Rb
end

Re_w, Rf_w, Rb_w = Rb_for(:water, 0)

println("=== Effect of fluid mixture on Rb — single U-loop, first order, T = $(T0)°C, Q = 15 L/min ===")
println(rpad("Fluid", 8), " ", rpad("% wt", 5), " ", rpad("Re", 8), " ", rpad("Rf [mK/W]", 10), " ",
        rpad("Rb [mK/W]", 10), " ", "ΔRb vs water")
println(rpad("water", 8), " ", rpad("-", 5), " ", rpad(round(Int, Re_w), 8), " ",
        rpad(round(Rf_w, digits=4), 10), " ", rpad(round(Rb_w, digits=4), 10), " ", "-")

for fluid in (:MPG, :MEG, :MMA, :MEA, :MKA, :MKF), percentage in (15, 30, 45)
    Re, Rf, Rb = Rb_for(fluid, percentage)
    ΔRb = 100 * (Rb - Rb_w) / Rb_w
    println(rpad(fluid, 8), " ", rpad(percentage, 5), " ", rpad(round(Int, Re), 8), " ",
            rpad(round(Rf, digits=4), 10), " ", rpad(round(Rb, digits=4), 10), " ",
            "+" * string(round(ΔRb, digits=1)) * "%")
end
