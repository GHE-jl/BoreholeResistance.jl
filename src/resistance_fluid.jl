"""
    Reynolds(V̇, r, ρf, μf)

Reynolds number of a fluid flowing in a pipe. Parameters `ρf` and `μf` depend on the fluid
temperature, and can be computed with functions `water_ρ` and `water_μ` respectively. This function
is valid for both cylinder pipes and annulus regions.
# Arguments
    - `V̇`: Fluid *speed* in pipe [m/s]
    - `r`: Pipe inside or annulus (r = rb - ro) radius [m]
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
# Output
    - `Re`: Reynolds number [-]
"""
function Reynolds(V̇::Real, r::Real, ρf::Real, μf::Real)
    return 2 * r * ρf * V̇ / μf
end

"""
    Prandtl(kf, cf, μf)
    Prandtl(kf, Cf, ρf, μf)

Prandtl number of a fluid flowing in a pipe. Parameters `kf`, `cf`, and `μf` depend on the fluid
temperature, and can be computed with functions `water_cp`, `water_μ` and `water_k` respectively.
# Arguments
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `Cf`: Fluid volumetric specific heat [J/m³K] (`water_cp(T) * water_ρ(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
# Output
    - `Pr`: Prandtl number [-]
"""
function Prandtl(kf::Real, cf::Real, μf::Real)
    return cf * μf / kf
end
function Prandtl(kf::Real, Cf::Real, ρf::Real, μf::Real)
    return Cf * μf / (kf * ρf)
end

"""
    friction_factor_Colebrook_White(Re, r, ϵ)

Darcy friction factor via the implicit Colebrook-White equation, solved iteratively. Valid for
laminar flow (Re < 2300, returns 64/Re) and turbulent flow. Valid for both cylinder pipes and
annulus regions.
# Arguments
    - `Re`: Reynolds number [-]
    - `r`: Pipe inside or annulus (r = rb - ro) radius [m]
    - `ϵ`: Pipe roughness [m]
# Output
    - `f`: Friction factor [-]
"""
function friction_factor_Colebrook_White(Re::Real, r::Real, ϵ::Real)
    if Re < eps()
        return 0.0
    elseif Re < 2300
        return 64 / Re
    else
        f = 0.02
        err = 1.0
        while err > 1e-5
            f_ = f
            f = (1 / (-2 * log10(ϵ / (3.7 * (2 * r)) + 2.51 / (Re * sqrt(f)))))^2
            err = abs(f - f_)
        end
        return f
    end
end

"""
    friction_factor_Tkachenko_Mileikovskyi(Re, r, ϵ)

Darcy friction factor via the explicit Tkachenko-Mileikovskyi approximation. Agrees with the
Colebrook-White equation within 1% for typical GHE flow conditions. Valid for laminar flow
(Re < 2320, returns 64/Re) and turbulent flow.
# Arguments
    - `Re`: Reynolds number [-]
    - `r`: Pipe inside or annulus (r = rb - ro) radius [m]
    - `ϵ`: Pipe roughness [m]
# Output
    - `f`: Friction factor [-]
# Reference
    - Mileikovskyi, V., & Tkachenko, T. (2021). Precise Explicit Approximations of the
        Colebrook-White Equation for Engineering Systems. In Z. Blikharskyy (Ed.), Proceedings of 
        EcoComfort 2020 (pp. 303–310). Springer International Publishing.
        https://doi.org/10.1007/978-3-030-57340-9_37
"""
function friction_factor_Tkachenko_Mileikovskyi(Re::Real, r::Real, ϵ::Real)
    if Re < eps()
        return 0.0
    elseif Re < 2320
        return 64 / Re
    end
    A₀ = -0.79638 * log(ϵ / (2 * r * 8.208) + 7.3357 / Re)
    A₁ = Re * ϵ / (2 * r) + A₀*9.3120665
    f = ((8.128943 + A₁) / (8.128943 * A₀ - 0.86859209 * A₁ * log(A₁ / (3.7099535 * Re))))^2
    return f
end

"""
    Nusselt(Re, Pr, r, ϵ=5e-6)
    Nusselt(V̇, r, kf, cf, ρf, μf, ϵ=5e-6)

Nusselt number of a fluid flowing in a pipe. The function is based on the Gnielinski correlation,
which is valid for laminar, transition and turbulent flow. The function assumes that the fluid is
flowing in a cylinder pipe, or in an annulus region.
# Arguments
    - `Re`: Reynolds number [-]
    - `Pr`: Prandtl number [-]
    - `V̇`: Fluid *speed* in pipe [m/s]
    - `r`: Pipe inside or annulus (r = rb - ro) radius [m]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 5e-6 for HDPE pipes)
# Output
    - `Nu`: Nusselt number [-]
"""
function Nusselt(Re::Real, Pr::Real, r::Real, ϵ::Real=5e-6)
    if Re < 2300                                # Laminar phase
        # Average between 4.364 for UHF and 3.657 for UBW boundary conditions
        return 4.0                              # Eq. 2.42 of Lamarche 2023
    elseif Re >= 2300 && Re < 4000              # Transition between laminar and turbulent flow
        γ = (Re - 2300) / (4000 - 2300)         # Eq. 2.49 of Lamarche 2023
        f = friction_factor_Colebrook_White(Re, r, ϵ)
        # Gnielinski (Eq. 2.43b of Lamarche 2023)
        Nu_4k = (f / 8) * (4000 - 1000) * Pr / (1 + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
        return (1 - γ) * 4 + γ * Nu_4k          # Eq. 2.48 of Lamarche 2023
    elseif Re >= 4000                           # Turbulent flow in a pipe
        f = friction_factor_Colebrook_White(Re, r, ϵ)
        # Gnielinski (Eq. 2.43b of Lamarche 2023)
        return (f / 8) * (Re - 1000) * Pr / (1 + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
    else
        error("Reynolds number must be non-negative.")
    end
end
function Nusselt(V̇::Real, r::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=5e-6)
    Re = Reynolds(V̇, r, ρf, μf)
    Pr = Prandtl(kf, cf, μf)
    return Nusselt(Re, Pr, r, ϵ)
end

"""
    Nusselt_annulus(Re, Pr, rb, ro, ϵo=5e-6, ϵi=5e-6)
    Nusselt_annulus(V̇, rb, ro, kf, cf, ρf, μf, ϵo=5e-6, ϵi=5e-6)

Function that computes the Nusselt number of a fluid flowing in an annulus region. The function is
based on the Gnielinski correlation, which is valid for laminar, transition and turbulent flow.
# Arguments
    - `Re`: Reynolds number [-]
    - `Pr`: Prandtl number [-]
    - `V̇`: Fluid *speed* in pipe [m/s]
    - `rb`: Outer radius of the annulus region [m]
    - `ro`: Inner radius of the annulus region [m]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
    - `ϵo`: Outer pipe roughness [m] (default 5e-6 for HDPE pipes)
    - `ϵi`: Inner pipe roughness [m] (default 5e-6 for HDPE pipes)
    # Output
    - `Nu`: Nusselt number [-]
# Reference
    - Lamarche, L. (2021). Analytic models and effective resistances for coaxial ground heat
        exchangers. Geothermics, 97, 102224. https://doi.org/10.1016/j.geothermics.2021.102224
"""
function Nusselt_annulus(Re::Real, Pr::Real, rb::Real, ro::Real, ϵo::Real=5e-6, ϵi::Real=5e-6)
    a = ro / rb
    if Re < 2300                                # Laminar phase
        return 3.66 + 1.2 * sqrt(a)             # Eq. 64b of Lamarche 2021
    elseif Re >= 2300 && Re < 4000              # Transition between laminar and turbulent flow
        Fₐ = (0.75 * a^-0.17 + (0.9 - 0.15 * a^0.6)) / (1 + a) # Eq. 64a of Lamarche 2021
        k₁ = 1.07 + (900 / 4000) - (0.63 / (1 + 10 * Pr))
        ϵ = (ϵo * rb + ϵi * ro) / (rb + ro)     # Equivalent roughness of the annulus region
        f = friction_factor_Colebrook_White(4000, rb - ro, ϵ)
        Nu_4k = Fₐ * (f / 8) * (4000 - 1000) * Pr / (k₁ + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
        γ = (Re - 2300) / (4000 - 2300)
        return (1 - γ) * (3.66 + 1.2 * sqrt(a)) + γ * Nu_4k
    elseif Re >= 4000
        Fₐ = (0.75 * a^-0.17 + (0.9 - 0.15 * a^0.6)) / (1 + a) # Eq. 64a of Lamarche 2021
        k₁ = 1.07 + (900 / Re) - (0.63 / (1 + 10 * Pr))
        ϵ = (ϵo * rb + ϵi * ro) / (rb + ro)     # Equivalent roughness of the annulus region
        f = friction_factor_Colebrook_White(Re, rb - ro, ϵ)
        return Fₐ * (f / 8) * (Re - 1000) * Pr / (k₁ + (12.7 * (f / 8)^0.5 * (Pr^(2 / 3) - 1)))
    else
        error("Reynolds number must be non-negative.")
    end
end
function Nusselt_annulus(V̇::Real, rb::Real, ro::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    ϵo::Real=5e-6, ϵi::Real=5e-6)
    Re = Reynolds(V̇, rb - ro, ρf, μf)
    Pr = Prandtl(kf, cf, μf)
    return Nusselt_annulus(Re, Pr, rb, ro, ϵo, ϵi)
end

"""
    convection_coefficient(Nu, r, kf)
    convection_coefficient(V̇, r, kf, cf, ρf, μf, ϵ=5e-6)

Convective heat transfer coefficient of a fluid flowing in a single cylinder pipe.
# Arguments
    - `Nu`: Nusselt number [-]
    - `V̇`: Fluid *speed* in pipe [m/s]
    - `r`: Pipe inner radius [m]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 5e-6 for HDPE pipes)
# Output
    - `h`: Convective heat transfer coefficient [W/m²K]
# Reference
    - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application.
        Springer Nature Switzerland.
"""
function convection_coefficient(Nu::Real, r::Real, kf::Real)
    return Nu * kf / (2 * r)                    # Eq. 2.32 of Lamarche 2023
end
function convection_coefficient(V̇::Real, r::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=5e-6)
    Nu = Nusselt(V̇, r, kf, cf, ρf, μf, ϵ)   # Eq. 2.42 of Lamarche 2023
    return convection_coefficient(Nu, r, kf)
end
"""
    resistance_fluid(Nu, r, kf)
    resistance_fluid(V̇, r, kf, cf, ρf, μf, ϵ=5e-6)

Convective thermal resistance of a fluid flowing in a single cylinder pipe.
# Arguments
    - `Nu`: Nusselt number [-]
    - `V̇`: Fluid *speed* in pipe [m/s]
    - `r`: Pipe inner radius [m]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid viscosity [kg/m⋅s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 5e-6 for HDPE pipes)
# Output
    - `Rf`: Fluid convective thermal resistance [mK/W]
# Reference
    - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application.
        Springer Nature Switzerland.
"""
function resistance_fluid(Nu::Real, r::Real, kf::Real)
    h = convection_coefficient(Nu, r, kf)
    return 1 / (2 * pi * r * h)                # Eq. 5.6 of Lamarche 2023
end
function resistance_fluid(V̇::Real, r::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=5e-6)
    Nu = Nusselt(V̇, r, kf, cf, ρf, μf, ϵ)
    return resistance_fluid(Nu, r, kf)
end
