using Test
using BoreholeResistance

@testset "BoreholeResistance.jl" begin

    @testset "Water properties" begin
        # Reference values at 20 °C (Engineering Toolbox)
        @test water_k(20.0)  ≈ 0.5984   rtol=0.005
        @test water_cp(20.0) ≈ 4182.0   rtol=0.005
        @test water_ρ(20.0)  ≈ 998.2    rtol=0.005
        @test water_μ(20.0)  ≈ 1.002e-3 rtol=0.01

        # Physical monotonicity over 10–90 °C
        Ts = 10.0:10.0:90.0
        @test all(diff([water_k(T) for T in Ts]) .> 0)   # conductivity rises with T
        @test all(diff([water_ρ(T) for T in Ts]) .< 0)   # density falls with T (above 4 °C)
        @test all(diff([water_μ(T) for T in Ts]) .< 0)   # viscosity falls with T
    end

    @testset "Reynolds number" begin
        # Re = 2r ρf V̇ / μf  — exact formula check
        @test Reynolds(1.0, 0.01, 1000.0, 0.001) ≈ 20000.0
        @test Reynolds(0.0, 0.01, 1000.0, 0.001) ≈ 0.0
    end

    @testset "Prandtl number" begin
        kf, cf, ρf, μf = 0.598, 4182.0, 998.2, 1.002e-3
        # Two overloads must be identical: cf·μf/kf == Cf·μf/(kf·ρf)
        @test Prandtl(kf, cf, μf) ≈ Prandtl(kf, cf * ρf, ρf, μf)
        # Water at 20 °C: Pr ≈ 7
        @test Prandtl(water_k(20.0), water_cp(20.0), water_μ(20.0)) ≈ 7.0 rtol=0.05
    end

    @testset "Nusselt number" begin
        r = 0.015
        # Laminar regime returns exactly 4.0 (hard-coded in Gnielinski correlation)
        @test Nusselt(1000.0, 7.0, r) == 4.0
        # Turbulent Nu must exceed the laminar value
        @test Nusselt(20000.0, 7.0, r) > 4.0
        # Two overloads must agree
        kf, cf, ρf, μf = water_k(20.0), water_cp(20.0), water_ρ(20.0), water_μ(20.0)
        V̇ = 1.0
        Re = Reynolds(V̇, r, ρf, μf)
        Pr = Prandtl(kf, cf, μf)
        @test Nusselt(V̇, r, kf, cf, ρf, μf) ≈ Nusselt(Re, Pr, r)
    end

    @testset "Friction factor" begin
        r, ϵ = 0.015, 1e-5
        # Laminar: f = 64/Re exactly
        @test friction_factor_Colebrook_White(1000.0, r, ϵ)        ≈ 64.0 / 1000.0
        @test friction_factor_Tkachenko_Mileikovskyi(1000.0, r, ϵ) ≈ 64.0 / 1000.0
        # Zero flow: f = 0
        @test friction_factor_Colebrook_White(0.0, r, ϵ)        ≈ 0.0
        @test friction_factor_Tkachenko_Mileikovskyi(0.0, r, ϵ) ≈ 0.0
        # Both methods must agree in turbulent regime (within 0.1 %)
        @test friction_factor_Colebrook_White(20000.0, r, ϵ) ≈
              friction_factor_Tkachenko_Mileikovskyi(20000.0, r, ϵ)  rtol=0.001
        # Turbulent f is positive and decreasing with Re
        f_turb = friction_factor_Colebrook_White(20000.0, r, ϵ)
        @test f_turb > 0
        @test friction_factor_Colebrook_White(40000.0, r, ϵ) < f_turb
    end

    @testset "Pipe resistance" begin
        ro, ri, kp = 0.020, 0.0164, 0.4
        Rp = resistance_pipe(ro, ri, kp)
        # Exact analytical formula
        @test Rp ≈ log(ro / ri) / (2π * kp)
        # Physical monotonicity
        @test resistance_pipe(ro, ri, 2 * kp) < Rp  # higher kp → lower Rp
        @test resistance_pipe(2 * ro, ri, kp)  > Rp  # thicker wall → higher Rp
    end

    @testset "Fluid resistance" begin
        r, kf = 0.0164, water_k(10.0)
        Nu = 4.0   # laminar
        # h = Nu·kf/(2r)  →  Rf = 1/(2πr·h) = 1/(π·Nu·kf)
        @test resistance_fluid(Nu, r, kf) ≈ 1.0 / (π * Nu * kf)
        # Higher Nu (better convection) → lower resistance
        @test resistance_fluid(20.0, r, kf) < resistance_fluid(4.0, r, kf)
        # Two overloads must agree
        kf, cf, ρf, μf = water_k(10.0), water_cp(10.0), water_ρ(10.0), water_μ(10.0)
        V̇ = 0.001
        Nu_num = Nusselt(V̇, r, kf, cf, ρf, μf)
        @test resistance_fluid(V̇, r, kf, cf, ρf, μf) ≈ resistance_fluid(Nu_num, r, kf)
    end

    # Parameters from script/script_single_Uloop.jl — typical HDPE single U-loop
    rb, ro, ri, s   = 0.075, 0.020, 0.0164, 0.08
    ks, kg, kp, ϵ   = 2.0, 1.0, 0.4, 5e-6
    H, T0           = 100.0, 10.0
    kf, cf, ρf, μf  = water_k(T0), water_cp(T0), water_ρ(T0), water_μ(T0)
    V_nom           = 15.0 / 1000 / 60   # 15 L/min → m³/s per pipe
    V̇_nom          = V_nom / (π * ri^2)  # mean fluid speed [m/s]
    Rp              = resistance_pipe(ro, ri, kp)
    Rf_nom          = resistance_fluid(V̇_nom, ri, kf, cf, ρf, μf, ϵ)

    @testset "Single U-loop — 15 L/min" begin
        @test Rp     > 0
        @test Rf_nom > 0

        for order in [0, 1]
            Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
            Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
            Rg = Rb - (Rp + Rf_nom) / 2          # grout resistance, N = 2 (Eq. 3)

            # Positivity
            @test Rb > 0
            @test Ra > 0
            @test Rg > 0

            # Both function signatures must give identical results
            @test Rb ≈ resistance_ULoop_borehole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; order=order)
            @test Ra ≈ resistance_ULoop_total_internal(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; order=order)
        end

        # Reference check against Javed and Spitler (2017), Table 6 (first-order column).
        # θ2 = 3 (rb = 3·rpo), moderate spacing θ1 = 0.445, ground λ = 3, grout 0.6, Rp = 0.05.
        rpo_t6 = 0.016; rb_t6 = 3 * rpo_t6; s_t6 = 2 * rb_t6 * 0.445
        Rb_t6 = resistance_ULoop_borehole(s_t6, rb_t6, rpo_t6, 3.0, 0.6, 0.05, 0.0; order=1)
        @test (Rb_t6 - 0.05 / 2) ≈ 0.14297 atol=2e-3     # Rg = Rb - Rp/N

        # First-order must converge toward zeroth-order (within 20 % for this geometry)
        Rb0 = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; order=0)
        Rb1 = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
        @test Rb0 ≈ Rb1 rtol=0.20

        # Effective borehole resistance — single U-loop uses V_nom as total system flow
        Ra1  = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
        Rbₑ  = resistance_ULoop_effective(V_nom, H, cf, ρf, Rb1, Ra1)
        @test Rbₑ >= Rb1   # short-circuiting can only increase the effective resistance

        # All three overloads of resistance_ULoop_effective must agree
        Rbₑ2 = resistance_ULoop_effective(V_nom, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_nom)
        Rbₑ3 = resistance_ULoop_effective(V_nom, H, s, rb, ro, ri, ks, kg, kp, kf,
            cf, ρf, μf, ϵ)
        @test Rbₑ ≈ Rbₑ2
        @test Rbₑ ≈ Rbₑ3

        # Longer borehole → larger temperature gradient along legs → more short-circuiting
        Rbₑ_short = resistance_ULoop_effective(V_nom, 50.0, cf, ρf, Rb1, Ra1)
        @test Rbₑ_short <= Rbₑ
    end

    @testset "Double U-loop — Claesson and Javed (2019) example" begin
        # Eq. 48 / Tables 1-3: λ = 3, λb = 1.5, rb = 57.5 mm, rp = 16 mm, H = 200 m, Rp = 0.05.
        # V is the flow in one U-tube loop (total system flow 2·Vf = 1.50 m³/h).
        ks_d, kg_d, rb_d, ro_d = 3.0, 1.5, 0.0575, 0.016
        H_d, Rp_d, cf_d, ρf_d  = 200.0, 0.05, 4180.0, 997.0
        V_loop = 1.50 / 3600 / 2

        # Moderate spacing, first-order column (rc = 28.11 mm → s = 2·rc):
        #   Rb, Ra_diag, Ra_adj, Rbe_dUHF, Rbe_dUBW, Rbe_aUHF, Rbe_aUBW
        rc = 0.02811; s_d = 2 * rc
        Rb_r, Rad_r, Raa_r = 0.06047, 0.20654, 0.31142
        dU, dW, aU, aW      = 0.10329, 0.09824, 0.08887, 0.08652

        Rb  = resistance_ULoop_borehole(s_d, rb_d, ro_d, ks_d, kg_d, Rp_d, 0.0; nLoop=2, order=1)
        Rad = resistance_ULoop_total_internal(s_d, rb_d, ro_d, ks_d, kg_d, Rp_d, 0.0; nLoop=2,
                  order=1, network="diagonal")
        Raa = resistance_ULoop_total_internal(s_d, rb_d, ro_d, ks_d, kg_d, Rp_d, 0.0; nLoop=2,
                  order=1, network="adjacent")

        # Borehole and internal resistance must match Table 2 (first-order column)
        @test Rb  ≈ Rb_r  atol=1e-3
        @test Rad ≈ Rad_r atol=1e-3
        @test Raa ≈ Raa_r atol=1e-3

        # Effective resistance for both flow networks and both boundary conditions (Tables/Eq. 44,46)
        @test resistance_ULoop_effective(V_loop, H_d, cf_d, ρf_d, Rb, Rad; nLoop=2, model="UHF") ≈ dU atol=1e-3
        @test resistance_ULoop_effective(V_loop, H_d, cf_d, ρf_d, Rb, Rad; nLoop=2, model="UBW") ≈ dW atol=1e-3
        @test resistance_ULoop_effective(V_loop, H_d, cf_d, ρf_d, Rb, Raa; nLoop=2, model="UHF") ≈ aU atol=1e-3
        @test resistance_ULoop_effective(V_loop, H_d, cf_d, ρf_d, Rb, Raa; nLoop=2, model="UBW") ≈ aW atol=1e-3

        # Diagonal internal resistance is smaller than adjacent (paper discussion)
        @test Rad < Raa

        # Core and geometry-based effective overloads must agree (nLoop = 2)
        Rbe_core = resistance_ULoop_effective(V_loop, H_d, cf_d, ρf_d, Rb, Rad; nLoop=2)
        Rbe_ovl  = resistance_ULoop_effective(V_loop, H_d, s_d, rb_d, ro_d, ks_d, kg_d, cf_d, ρf_d,
                       Rp_d, 0.0; nLoop=2, network="diagonal")
        @test Rbe_core ≈ Rbe_ovl
        @test Rbe_core >= Rb   # short-circuiting can only increase the effective resistance

        # Longer borehole → more short-circuiting → larger effective resistance
        @test resistance_ULoop_effective(V_loop, 50.0, cf_d, ρf_d, Rb, Rad; nLoop=2) <= Rbe_core
    end

    @testset "Double U-loop — HDPE overloads and single/double comparison" begin
        for order in [0, 1]
            Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=order)
            Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2,
                order=order)
            @test Rb > 0
            @test Ra > 0
            # Long- and short-form overloads must agree (long form receives per-loop flow V_nom)
            @test Rb ≈ resistance_ULoop_borehole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; nLoop=2, order=order)
            @test Ra ≈ resistance_ULoop_total_internal(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; nLoop=2, order=order)
        end

        # At equal per-pipe flow (same Rf) the double U-loop has the lower borehole resistance
        Rb_single = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=1, order=1)
        Rb_double = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=1)
        @test Rb_double < Rb_single
    end

    # Coaxial exchanger — "Borehole 1" of Lamarche (2021), Table 1 (diameters → radii)
    @testset "Coaxial GHE — Lamarche 2021, Borehole 1" begin
        Hc  = 150.0
        rbc = 0.1524 / 2
        roo = 0.1143 / 2
        roi = 0.0973 / 2
        rio = 0.0603 / 2
        rii = 0.0494 / 2
        kpi, kpo, kgc = 0.4, 0.4, 1.7
        ϵc  = 5e-6
        T0c = 10.0
        kfc, cfc, ρfc, μfc = water_k(T0c), water_cp(T0c), water_ρ(T0c), water_μ(T0c)
        ṁc  = 3.25
        Vc  = ṁc / ρfc

        # Published resistances from the paper's convection coefficients (Section 2.3):
        # R'12 = 0.084 and R'1 = 0.092 mK/W (hin = 5170, hann = 2150 W/m²K)
        hin_paper, hann_paper = 5170.0, 2150.0
        R1p, R12p = resistance_coaxial(rii, rio, roi, roo, rbc, kgc, kpi, kpo, hin_paper,
            hann_paper)
        @test R12p ≈ 0.084 atol = 0.002   # Eq. 1 of Lamarche 2021
        @test R1p  ≈ 0.092 atol = 0.002   # Eq. 2 of Lamarche 2021

        # Effective resistances match Table 2 (Rb* ≈ 0.0915 UHF, 0.0916 UBW). Section 4 of the
        # paper neglects the convection resistances (h → ∞), giving R'1 = 0.091, R'12 = 0.0798.
        R1n, R12n = resistance_coaxial(rii, rio, roi, roo, rbc, kgc, kpi, kpo, 1e12, 1e12)
        @test R1n  ≈ 0.091  atol = 0.002
        @test R12n ≈ 0.0798 atol = 0.002
        Rbe_UHF = resistance_coaxial_effective(Vc, Hc, cfc, ρfc, R1n, R12n; model="UHF")
        Rbe_UBW = resistance_coaxial_effective(Vc, Hc, cfc, ρfc, R1n, R12n; model="UBW")
        @test Rbe_UHF ≈ 0.0915 atol = 0.001   # Eq. 31
        @test Rbe_UBW ≈ 0.0916 atol = 0.001   # Eq. 14
        # Mean model lies between the two
        Rbe_mean = resistance_coaxial_effective(Vc, Hc, cfc, ρfc, R1n, R12n; model="mean")
        @test min(Rbe_UHF, Rbe_UBW) <= Rbe_mean <= max(Rbe_UHF, Rbe_UBW)

        # Short-circuiting can only increase the effective resistance beyond Rb = R1
        @test Rbe_UHF >= R1n
        @test Rbe_UBW >= R1n

        # Built-in Nusselt correlations give positive, physically sensible resistances
        R1, R12 = resistance_coaxial(Vc, rii, rio, roi, roo, rbc, kgc, kpi, kpo, kfc, cfc, ρfc,
            μfc, ϵc)
        @test R1  > 0
        @test R12 > 0

        # Both overloads of the effective resistance must agree
        Rbe_short = resistance_coaxial_effective(Vc, Hc, cfc, ρfc, R1, R12; model="UHF")
        Rbe_full  = resistance_coaxial_effective(Vc, Hc, rii, rio, roi, roo, rbc, kgc,
            kpi, kpo, kfc, cfc, ρfc, μfc, ϵc; model="UHF")
        @test Rbe_short ≈ Rbe_full

        # Longer borehole / lower flow → more short-circuiting → larger Rb*
        Rbe_long  = resistance_coaxial_effective(Vc, 300.0, cfc, ρfc, R1n, R12n)
        @test Rbe_long >= Rbe_UHF
        Rbe_slow  = resistance_coaxial_effective(Vc / 4, Hc, cfc, ρfc, R1n, R12n)
        @test Rbe_slow >= Rbe_UHF

        # Unknown model must error
        @test_throws ErrorException resistance_coaxial_effective(Vc, Hc, cfc, ρfc, R1n,
            R12n; model="foo")
    end

    # "UHF_gradient" model (Eqs. 58-63) — Borehole 3b of Lamarche (2021), Table 3
    @testset "Coaxial GHE — Lamarche 2021, Borehole 3b, UHF_gradient" begin
        Hg  = 200.0
        R1g, R12g = 0.0209, 0.0798        # Table 2 (Borehole 3b) and Section 4 (Borehole 3)
        T0g = 10.0
        cfg, ρfg = water_cp(T0g), water_ρ(T0g)
        ṁg  = 0.8
        Vg  = ṁg / ρfg

        Rbe_UHF  = resistance_coaxial_effective(Vg, Hg, cfg, ρfg, R1g, R12g; model="UHF")
        Rbe_grad = resistance_coaxial_effective(Vg, Hg, cfg, ρfg, R1g, R12g; model="UHF_gradient")
        @test Rbe_UHF  ≈ 0.0356 atol = 0.001   # Eq. 31, gradient-independent by construction
        @test Rbe_grad ≈ 0.0420 atol = 0.001   # Eq. 61/63, Table 3 case 1 (heat injection, annulus-in)

        # The gradient correction amplifies the short-circuit beyond the plain UHF estimate
        @test Rbe_grad > Rbe_UHF > R1g

        # It is a distinct closed form, not a generalization of UHF (does not collapse onto it)
        @test !isapprox(Rbe_grad, Rbe_UHF; rtol=0.05)

        # Full-geometry overload must agree with the short-form overload (same R1, R12)
        rbc3, roo3, roi3, rio3, rii3 = 0.1524 / 2, 0.1143 / 2, 0.0973 / 2, 0.0603 / 2, 0.0494 / 2
        R1full, R12full = resistance_coaxial(Vg, rii3, rio3, roi3, roo3, rbc3, 2.5, 0.4, 10.0,
            water_k(T0g), cfg, ρfg, water_μ(T0g), 5e-6)
        Rbe_grad_short = resistance_coaxial_effective(Vg, Hg, cfg, ρfg, R1full, R12full;
            model="UHF_gradient")
        Rbe_grad_full  = resistance_coaxial_effective(Vg, Hg, rii3, rio3, roi3, roo3, rbc3, 2.5,
            0.4, 10.0, water_k(T0g), cfg, ρfg, water_μ(T0g), 5e-6; model="UHF_gradient")
        @test Rbe_grad_short ≈ Rbe_grad_full
    end

    # testing the density calculation of fluid_property
    @testset "fluid density (fluid_property) - Enginering Toolbox" begin

        # for 10% propylene glycol
        ρ1_ref = 1012 # 0 °C
        ρ2_ref = 998 # 40 °C
        ρ3_ref = 976 # 80 °C
        ρ4_ref = 965 # 100 °C

        # for 30% propylene glycol
        ρ5_ref = 1031 # 0 °C
        ρ6_ref = 1010 # 40 °C
        ρ7_ref = 983 # 80 °C
        ρ8_ref = 969 # 100 °C

        # for 60% propylene glycol
        ρ9_ref = 1061 # 0 °C
        ρ10_ref = 1029 # 40 °C
        ρ11_ref = 994 # 80 °C
        ρ12_ref = 976 # 100 °C


        _,_, ρ1, _ = fluid_property(0.0, :MPG; percentage= 10)
        _,_, ρ2, _ = fluid_property(40.0, :MPG; percentage= 10)
        _,_, ρ3, _ = fluid_property(80.0, :MPG; percentage= 10)
        _,_, ρ4, _ = fluid_property(100.0, :MPG; percentage= 10)

        _,_, ρ5, _ = fluid_property(0.0, :MPG; percentage= 30)
        _,_, ρ6, _ = fluid_property(40.0, :MPG; percentage= 30)
        _,_, ρ7, _ = fluid_property(80.0, :MPG; percentage= 30)
        _,_, ρ8, _ = fluid_property(100.0, :MPG; percentage= 30)

        _,_, ρ9, _ = fluid_property(0.0, :MPG; percentage= 60)
        _,_, ρ10, _ = fluid_property(40.0, :MPG; percentage= 60)
        _,_, ρ11, _ = fluid_property(80.0, :MPG; percentage= 60)
        _,_, ρ12, _ = fluid_property(100.0, :MPG; percentage= 60)


        @test ρ1_ref ≈ ρ1 rtol = 1e-2
        @test ρ2_ref ≈ ρ2 rtol = 1e-2
        @test ρ3_ref ≈ ρ3 rtol = 1e-2
        @test ρ4_ref ≈ ρ4 rtol = 1e-2

        @test ρ5_ref ≈ ρ5 rtol = 1e-2
        @test ρ6_ref ≈ ρ6 rtol = 1e-2
        @test ρ7_ref ≈ ρ7 rtol = 1e-2
        @test ρ8_ref ≈ ρ8 rtol = 1e-2

        @test ρ9_ref ≈ ρ9 rtol = 1e-2
        @test ρ10_ref ≈ ρ10 rtol = 1e-2
        @test ρ11_ref ≈ ρ11 rtol = 1e-2
        @test ρ12_ref ≈ ρ12 rtol = 1e-2
    end

    @testset "fluid_property(:water) vs water_ series" begin
        # Stops short of 100 °C: at P = 100 kPa water's actual boiling point is ≈99.6 °C, above
        # which CoolProp's equation of state switches to the vapor branch (see fluid_property.jl).
        for T in 0.0:10.0:90.0
            k, cp, ρ, μ = fluid_property(T, :water)
            @test k  ≈ water_k(T)  rtol = 5e-3
            @test cp ≈ water_cp(T) rtol = 5e-3
            @test ρ  ≈ water_ρ(T)  rtol = 5e-3
            @test μ  ≈ water_μ(T)  rtol = 2e-2
        end
    end

    @testset "fluid_property error handling" begin
        # Unsupported fluid symbol must error, never silently fall back to water
        @test_throws ErrorException fluid_property(20.0, :notafluid)

        # Out-of-composition-range mixture request must error (CoolProp's own bound)
        @test_throws Exception fluid_property(20.0, :MPG; percentage=90)

        # Fractional percentages are supported (not restricted to integers)
        k, cp, ρ, μ = fluid_property(20.0, :MPG; percentage=25.5)
        @test ρ > 0

        # Out-of-range water temperature warns but still returns CoolProp's value
        @test (@test_logs (:warn,) match_mode=:any fluid_property(150.0, :water)) !== nothing
    end

    @testset "fluid_property additional mixtures (sanity)" begin
        for f in (:MPG, :MEG, :MMA, :MEA, :MKA, :MKF)
            k, cp, ρ, μ = fluid_property(20.0, f; percentage=30)
            @test k > 0
            @test cp > 0
            @test ρ > 0
            @test μ > 0
        end
    end

end
