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
            Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
            Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=order)
            Rg = Rb - Rp - Rf_nom

            # Positivity
            @test Rb > 0
            @test Ra > 0
            @test Rg > 0

            # Both function signatures must give identical results
            @test Rb ≈ resistance_borehole_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; order=order)
            @test Ra ≈ resistance_total_internal_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; order=order)
        end

        # First-order must converge toward zeroth-order (within 20 % for this geometry)
        Rb0 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=0)
        Rb1 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
        @test Rb0 ≈ Rb1 rtol=0.20

        # Effective borehole resistance — single U-loop uses V_nom as total system flow
        Ra1  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; order=1)
        Rbₑ  = resistance_borehole_effective(V_nom, H, cf, ρf, Rb1, Ra1)
        @test Rbₑ >= Rb1   # short-circuiting can only increase the effective resistance

        # All three overloads of resistance_borehole_effective must agree
        Rbₑ2 = resistance_borehole_effective(V_nom, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_nom)
        Rbₑ3 = resistance_borehole_effective(V_nom, H, s, rb, ro, ri, ks, kg, kp, kf,
            cf, ρf, μf, ϵ)
        @test Rbₑ ≈ Rbₑ2
        @test Rbₑ ≈ Rbₑ3

        # Longer borehole → larger temperature gradient along legs → more short-circuiting
        Rbₑ_short = resistance_borehole_effective(V_nom, 50.0, cf, ρf, Rb1, Ra1)
        @test Rbₑ_short <= Rbₑ
    end

    @testset "Double U-loop — 30 L/min total" begin
        V_total = 2 * V_nom   # 30 L/min total system flow

        for order in [0, 1]
            Rb = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=order)
            Ra = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2,
                order=order)

            # Positivity
            @test Rb > 0
            @test Ra > 0

            # Both function signatures must give identical results.
            # The long-form function receives V_nom (per-pipe flow) to compute the same Rf.
            @test Rb ≈ resistance_borehole_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; nLoop=2, order=order)
            @test Ra ≈ resistance_total_internal_multipole(V_nom, s, rb, ro, ri, ks, kg, kp, kf,
                cf, ρf, μf, ϵ; nLoop=2, order=order)
        end

        # First-order must converge toward zeroth-order
        Rb0 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=0)
        Rb1 = resistance_borehole_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=1)
        @test Rb0 ≈ Rb1 rtol=0.20

        # Both diagonal and adjacent networks must give positive results
        Ra_diag = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2,
            order=1, network="diagonal")
        Ra_adj  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2,
            order=1, network="adjacent")
        @test Ra_diag > 0
        @test Ra_adj  > 0

        # Effective resistance uses the total system flow (doubled) to get the correct
        # thermal capacity flow rate for the two-loop circuit
        Ra1  = resistance_total_internal_multipole(s, rb, ro, ks, kg, Rp, Rf_nom; nLoop=2, order=1)
        Rbₑ  = resistance_borehole_effective(V_total, H, cf, ρf, Rb1, Ra1)
        @test Rbₑ >= Rb1

        # Medium overload (Rp and Rf provided, V is total flow) must agree with short-form
        Rbₑ2 = resistance_borehole_effective(V_total, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf_nom;
            nLoop=2)
        @test Rbₑ ≈ Rbₑ2

        # Longer borehole → larger temperature gradient along legs → more short-circuiting
        Rbₑ_short = resistance_borehole_effective(V_total, 50.0, cf, ρf, Rb1, Ra1)
        @test Rbₑ_short <= Rbₑ
    end

end
