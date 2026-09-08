# Validation of the double U-loop thermal-resistance functions against Claesson and Javed (2019),
# "Explicit multipole formulas and thermal network models for calculating thermal resistances of
# double U-pipe borehole heat exchangers", Science and Technology for the Built Environment 25(8),
# 980-992.
#
# The example of Eq. 48 (Figure 11, Tables 1-3) is reproduced: borehole resistance Rb (Eqs. 13-14),
# internal resistance for diagonal (Eqs. 18-19) and adjacent (Eqs. 22-23) inlet pipes, and the
# effective resistance for uniform heat flux (Eq. 44) and uniform wall temperature (Eq. 46).
#
# Geometry convention: the four pipes sit on a circle of radius rc (the pipe-centre radius). This
# package parameterises the double U-loop through the shank spacing `s`, with rc = s/2, so each
# reference case is run with s = 2·rc.

using BoreholeResistance

# Parameters — Eq. 48 of Claesson and Javed (2019)
ks = 3.0                    # Ground thermal conductivity λ  [W/mK]
kg = 1.5                    # Grout thermal conductivity λb  [W/mK]
rb = 0.05750                # Borehole radius [m]
ro = 0.016                  # Pipe outer radius rp [m]
H  = 200.0                  # Borehole length [m]
Rp = 0.05                   # Total fluid-to-pipe resistance for one pipe [mK/W]
cf = 4180.0                 # Fluid specific heat [J/kgK]
ρf = 997.0                  # Fluid density [kg/m³]
V_loop = 1.50 / 3600 / 2    # Flow in ONE U-tube loop [m³/s]  (total 2·Vf = 1.50 m³/h)

# Three spacings of Figure 11, given as pipe-centre radius rc [m]
cases = [("Close", 0.02263), ("Moderate", 0.02811), ("Wide", 0.04150)]

# Reference values from Tables 1-3, zeroth and first-order multipole results, keyed by spacing name.
# Rb, Ra_diagonal, Ra_adjacent, Rbe_diag_UHF, Rbe_diag_UBW, Rbe_adj_UHF, Rbe_adj_UBW
ref0 = Dict(
    "Close"    => (0.08386, 0.17017, 0.29857, 0.13583, 0.13035, 0.11348, 0.11158),
    "Moderate" => (0.06658, 0.21150, 0.33220, 0.10840, 0.10394, 0.09320, 0.09129),
    "Wide"     => (0.03493, 0.26288, 0.36763, 0.06857, 0.06348, 0.05899, 0.05621),
)
ref1 = Dict(
    "Close"    => (0.07509, 0.16039, 0.26166, 0.13023, 0.12353, 0.10889, 0.10619),
    "Moderate" => (0.06047, 0.20654, 0.31142, 0.10329, 0.09824, 0.08887, 0.08652),
    "Wide"     => (0.03143, 0.26273, 0.36144, 0.06509, 0.05955, 0.05590, 0.05278),
)

println("=== Double U-loop — Zeroth-order borehole & internal resistance vs Claesson and Javed (2019) ===")
println(rpad("spacing", 9), " ", rpad("Rb", 8), " ", rpad("Rb ref", 8), " ",
        rpad("Ra_diag", 8), " ", rpad("Ra_d ref", 9), " ", rpad("Ra_adj", 8), " ", "Ra_a ref")
for (name, rc) in cases
    local s = 2 * rc                        # rc = s/2
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=0)
    Rad = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=0,
              network="diagonal")
    Raa = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=0,
              network="adjacent")
    Rb_ref, Rad_ref, Raa_ref = ref0[name][1], ref0[name][2], ref0[name][3]

    @assert isapprox(Rb,  Rb_ref;  atol=1e-3)  "Rb must match Table ($name)"
    @assert isapprox(Rad, Rad_ref; atol=1e-3)  "Ra diagonal must match Table ($name)"
    @assert isapprox(Raa, Raa_ref; atol=1e-3)  "Ra adjacent must match Table ($name)"

    println(rpad(name, 9), " ", rpad(round(Rb, digits=5), 8), " ", rpad(Rb_ref, 8), " ",
            rpad(round(Rad, digits=5), 8), " ", rpad(Rad_ref, 9), " ",
            rpad(round(Raa, digits=5), 8), " ", Raa_ref)
end
println()

println("=== Double U-loop — First-order borehole & internal resistance vs Claesson and Javed (2019) ===")
println(rpad("spacing", 9), " ", rpad("Rb", 8), " ", rpad("Rb ref", 8), " ",
        rpad("Ra_diag", 8), " ", rpad("Ra_d ref", 9), " ", rpad("Ra_adj", 8), " ", "Ra_a ref")
for (name, rc) in cases
    local s = 2 * rc                        # rc = s/2
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1)
    Rad = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1,
              network="diagonal")
    Raa = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1,
              network="adjacent")
    Rb_ref, Rad_ref, Raa_ref = ref1[name][1], ref1[name][2], ref1[name][3]

    @assert isapprox(Rb,  Rb_ref;  atol=1e-3)  "Rb must match Table ($name)"
    @assert isapprox(Rad, Rad_ref; atol=1e-3)  "Ra diagonal must match Table ($name)"
    @assert isapprox(Raa, Raa_ref; atol=1e-3)  "Ra adjacent must match Table ($name)"

    println(rpad(name, 9), " ", rpad(round(Rb, digits=5), 8), " ", rpad(Rb_ref, 8), " ",
            rpad(round(Rad, digits=5), 8), " ", rpad(Rad_ref, 9), " ",
            rpad(round(Raa, digits=5), 8), " ", Raa_ref)
end
println()

println("=== Double U-loop — Zeroth-order effective resistance vs Claesson and Javed (2019) ===")
println(rpad("spacing", 9), " ", rpad("network", 9), " ", rpad("UHF", 8), " ", rpad("UHF ref", 8),
        " ", rpad("UBW", 8), " ", "UBW ref")
for (name, rc) in cases
    local s = 2 * rc
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=0)
    _, _, _, ref_dU, ref_dW, ref_aU, ref_aW = ref0[name]

    for (net, refs) in [("diagonal", (ref_dU, ref_dW)), ("adjacent", (ref_aU, ref_aW))]
        Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=0,
                 network=net)
        RbeU = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Ra; nLoop=2, model="UHF")
        RbeW = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Ra; nLoop=2, model="UBW")

        @assert isapprox(RbeU, refs[1]; atol=1e-3)  "Rbe UHF must match Table ($name, $net)"
        @assert isapprox(RbeW, refs[2]; atol=1e-3)  "Rbe UBW must match Table ($name, $net)"
        @assert RbeU >= Rb  "short-circuiting can only increase the effective resistance"

        println(rpad(name, 9), " ", rpad(net, 9), " ", rpad(round(RbeU, digits=5), 8), " ",
                rpad(refs[1], 8), " ", rpad(round(RbeW, digits=5), 8), " ", refs[2])
    end
end
println()

println("=== Double U-loop — First-order effective resistance vs Claesson and Javed (2019) ===")
println(rpad("spacing", 9), " ", rpad("network", 9), " ", rpad("UHF", 8), " ", rpad("UHF ref", 8),
        " ", rpad("UBW", 8), " ", "UBW ref")
for (name, rc) in cases
    local s = 2 * rc
    Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1)
    _, _, _, ref_dU, ref_dW, ref_aU, ref_aW = ref1[name]

    for (net, refs) in [("diagonal", (ref_dU, ref_dW)), ("adjacent", (ref_aU, ref_aW))]
        Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1,
                 network=net)
        RbeU = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Ra; nLoop=2, model="UHF")
        RbeW = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Ra; nLoop=2, model="UBW")

        @assert isapprox(RbeU, refs[1]; atol=1e-3)  "Rbe UHF must match Table ($name, $net)"
        @assert isapprox(RbeW, refs[2]; atol=1e-3)  "Rbe UBW must match Table ($name, $net)"
        @assert RbeU >= Rb  "short-circuiting can only increase the effective resistance"

        println(rpad(name, 9), " ", rpad(net, 9), " ", rpad(round(RbeU, digits=5), 8), " ",
                rpad(refs[1], 8), " ", rpad(round(RbeW, digits=5), 8), " ", refs[2])
    end
end
println()

# The paper notes the internal (and effective) resistance is larger for diagonal than for
# adjacent inlet pipes. The next blocs test this.
s = 2 * 0.02811
Rb  = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1)
Rad = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1,
    network="diagonal")
Raa = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, 0.0; nLoop=2, order=1,
    network="adjacent")
Rbe_d = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Rad; nLoop=2, model="mean")
Rbe_a = resistance_ULoop_effective(V_loop, H, cf, ρf, Rb, Raa; nLoop=2, model="mean")
println("=== Double U-loop — Diagonal vs adjacent internal resistance (order 1) ===")
println("  Ra_diag = $(round(Rad, digits=4)) mK/W, Ra_adj = $(round(Raa, digits=4)) mK/W")
println("  Rbe_diag = $(round(Rbe_d, digits=4)) mK/W,  Rbe_adj = $(round(Rbe_a, digits=4)) mK/W")
@assert Rad < Raa      "diagonal internal resistance is smaller than adjacent"
@assert Rbe_d > Rbe_a  "diagonal inlet pipes give the larger effective resistance"
println()

# HDPE example with water properties — single vs double comparison and overload cross-check
rb2, ro2, ri2 = 0.075, 0.020, 0.0164
s2, ks2, kg2  = 0.08, 2.0, 1.0
kp2, ϵ2, T0   = 0.4, 5e-6, 10.0
kf, cf2, ρf2, μf = fluid_property(T0, :water)
V2 = 15.0 / 1000 / 60                         # per-loop flow [m³/s]
Rp2 = resistance_pipe(ro2, ri2, kp2)
Rf2 = resistance_fluid(V2 / (π * ri2^2), ri2, kf, cf2, ρf2, μf, ϵ2)

# Long- and short-form overloads must agree
for order in [0, 1]
    Rb_s = resistance_ULoop_borehole(s2, rb2, ro2, ks2, kg2, Rp2, Rf2; nLoop=2, order=order)
    Rb_l = resistance_ULoop_borehole(V2, s2, rb2, ro2, ri2, ks2, kg2, kp2, kf, cf2, ρf2, μf, ϵ2;
        nLoop=2, order=order)
    @assert Rb_s ≈ Rb_l  "resistance_ULoop_borehole overloads must agree (order $order)"
end

Rb1 = resistance_ULoop_borehole(s2, rb2, ro2, ks2, kg2, Rp2, Rf2; nLoop=1, order=1)
Rb2 = resistance_ULoop_borehole(s2, rb2, ro2, ks2, kg2, Rp2, Rf2; nLoop=2, order=1)
println("=== Single vs double U-loop (HDPE, same per-loop flow, order 1) ===")
println("  Rb single = $(round(Rb1, digits=4)) mK/W")
println("  Rb double = $(round(Rb2, digits=4)) mK/W")
@assert Rb2 < Rb1  "at equal per-pipe flow the double U-loop has the lower borehole resistance"
