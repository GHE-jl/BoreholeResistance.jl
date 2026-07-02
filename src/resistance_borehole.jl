
"""
    resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1)
    resistance_ULoop_borehole(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0,
        nLoop=1, order=1)

Function that computes the thermal borehole resistance of a ground heat exchanger based on the
zeroth- or first-order multipole formula for a single U-loop of Hellström (1991) for grout thermal 
resistance. It is the equivalent of the line source theory. See Eq. 8.36 of the Hellström (1991).
Note: To obtain only the grout thermal resistance, set `Rp` and `Rf` as 0.0.
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `Rp`: Pipe thermal resistance [mK/W]
    - `Rf`: Fluid thermal resistance [mK/W]
    - `order`: Order of the multipole method (default 1)
    - `nLoop`: Number of loops in the ground heat exchanger (2 pipes per loop) (default 1)
# Output
    - `Rb`: Borehole thermal resistance [mK/W]
# Reference
    - Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    - Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    - Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for
        calculating thermal resistances of double U-pipe borehole heat exchangers. Science and
        Technology for the Built Environment, 25(8), 980–992.
        https://doi.org/10.1080/23744731.2019.1620565
"""
function resistance_ULoop_borehole(s::Real, rb::Real, ro::Real, ks::Real, kg::Real,
    Rp::Real, Rf::Real; nLoop::Int=1, order::Int=1)
    # Initial parameters
    Rₚ = Rp + Rf                                    # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ                             # β parameter for the multipole method
    σ = (kg - ks) / (kg + ks)                       # σ parameter for the multipole method

    if nLoop == 1                                   # Single U-loop
        # Compute θ₁ to θ₃
        θ₁ = s / (2 * rb)   # equal to D/rb
        θ₂ = rb / ro
            
        if order == 0                               # Zeroth-order multipole method (line-source)
            # Compute Rb with Eq. 12 from Javed and Spitler 2017
            Rb = (1 / (4 * π * kg)) * (β + log(θ₂ / (2 * θ₁ * (1 - θ₁^4)^σ)))
        elseif order == 1                           # First-order multipole method
            θ₃ = ro / s
            b₁ = (1 + β) / (1 - β)
            # Compute Rb with Eq. 13 from Javed and Spitler 2017
            tmp = (1 - θ₁^4)
            a = log(θ₂ / (2 * θ₁ * tmp^σ))
            b = θ₃^2 * (1 - (4 * σ * θ₁^4 / tmp))^2
            c = b₁ + (θ₃^2 * (1 + (16 * σ * θ₁^4 / (tmp^2))))
            Rb = (1 / (4 * π * kg)) * (β + a - (b / c))
        else
            error("Only order 0 and 1 are implemented for the multipole method.")
        end
    elseif nLoop == 2                               # Double U-loop
        if order == 0                               # Zeroth-order multipole method (line-source)
            Rb = Rₚ / 4 + (1 / (4 * π * kg)) * (log(rb^4 / (4 * ro * (s / 2)^3)) + σ * log(rb^8 / 
                (rb^8 - (s / 2)^8)))
        elseif order == 1                           # First-order multipole method
            θ₁ = ro^2 / (4 * (s / 2)^2)
            θ₂ = (s / 2)^2 / (rb^8 - (s / 2)^8)^(1/4)
            θ₃ = rb^2 / (rb^8 - (s / 2)^8)^(1/4)
            b₁ = (1 - β) / (1 + β)

            Rb = Rₚ / 4 + (1 / (4 * π * kg)) * (log(rb^4 / (4 * ro * (s / 2)^3)) + σ * log(rb^8 / 
                (rb^8 - (s / 2)^8))) - ((b₁ * θ₁ * (3 - 8 * σ * θ₂^4)^2) / ((8 * π * kg) * 
                (1 + b₁ * θ₁ * (5 + 64 * σ * θ₂^4 * θ₃^4))))
        else
            error("Only order 0 and 1 are implemented for the multipole method.")
        end
    else
        error("Only single and double U-loops are implemented for the multipole method.")
    end

    return Rb
end
function resistance_ULoop_borehole(V::Real, s::Real, rb::Real, ro::Real, ri::Real, ks::Real,
    kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0;
    nLoop::Int=1, order::Int=1)
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp)

    # Compute Rb with the first-order multipole method
    return resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop, order=order)
end

"""
    resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf, nLoop=1, order=1,
        network="diagonal")
    resistance_ULoop_total_internal(V, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, ϵ=0.0,
        nLoop=1, order=1; network="diagonal")

Computes the zeroth- or first-order multipole method for the total internal resistance (Eq. 25 or 26
of Javed and Spitler 2017) by Hellström 1991 for a single U-tube ground heat exchanger. This is used
to convert the borehole thermal resistance (Rb) in effective borehole thermal resistance (named Rb*
or Rbₑ).
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `Rp`: Pipe thermal resistance [mK/W]
    - `Rf`: Fluid thermal resistance [mK/W]
    - `nLoop`: Number of loops in the ground heat exchanger (2 pipes per loop) (default 1)
    - `order`: Order of the multipole method (default 1)
    - `network`: `diagonal` and `adjacent` network for double U-loop (default `diagonal`)
# Output
    - `Ra`: Total internal thermal resistance [mK/W]
# Reference
    - Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods
        for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. 
        https://doi.org/10.1016/j.apenergy.2016.11.079
    - Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.”
        http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
    - Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for
        calculating thermal resistances of double U-pipe borehole heat exchangers. Science and
        Technology for the Built Environment, 25(8), 980–992.
        https://doi.org/10.1080/23744731.2019.1620565
"""
function resistance_ULoop_total_internal(s::Real, rb::Real, ro::Real, ks::Real, kg::Real,
    Rp::Real, Rf::Real; nLoop::Int=1, order::Int=1, network::String="diagonal")
    # Initial parameters
    Rₚ = Rp + Rf                                    # Pipe and fluid resistance
    β = 2 * π * kg * Rₚ                             # β parameter for the multipole method
    σ = (kg - ks) / (kg + ks)                       # σ parameter for the multipole method

    if nLoop == 1                                   # Single U-loop
        # Compute θ₁ and θ₃
        θ₁ = s / (2 * rb)
        θ₃ = ro / s

        # Eq. 26 from Javed and Spitler 2017
        if order == 0
            Ra = (1 / (π * kg)) * (β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ)))
        elseif order == 1
            Ra = (1 / (π * kg)) * (β + log((1 + θ₁^2)^σ / (θ₃ * (1 - θ₁^2)^σ)) -
                ((θ₃^2 * (1 - (θ₁^4) + (4 * σ * θ₁^2))^2) /
                (((1 + β) / (1 - β)) * (1 - θ₁^4)^2 - (θ₃^2 * (1 - θ₁^4)^2) +
                (8 * σ * θ₁^2 * θ₃^2 * (1 + θ₁^4)))))
        else
            error("Only order 0 and 1 are implemented for the multipole method.")
        end
    elseif nLoop == 2                               # Double U-loop
        if network == "diagonal"
            if order == 0
                # Eq. 18 from Claesson and Javed 2019 (diagonal network, zeroth-order multipole)
                Ra = 2 * Rₚ + (1 / (π * kg)) * (log(s / (2 * ro)) + σ * log((rb^4 + (s / 2)^4) / 
                    (rb^4 - (s / 2)^4)))
            elseif order == 1
                # Eq. 19 from Claesson and Javed 2019 (diagonal network, first-order multipole)
                θ₁ = ro^2 / (4 * (s / 2)^2)
                θ₂ = (s / 2)^2 / (rb^8 - (s / 2)^8)^(1/4)
                θ₃ = rb^2 / (rb^8 - (s / 2)^8)^(1/4)
                b₁ = (1 - β) / (1 + β)

                Ra = 2 * Rₚ + (1 / (π * kg)) * (log(s / (2 * ro)) + σ * log((rb^4 + (s / 2)^4) /
                    (rb^4 - (s / 2)^4))) - ((1 / (π * kg)) * (b₁ * θ₁ * (1 + 8 * σ * θ₂^2 * 
                    θ₃^2)^2) / (1 - b₁ * θ₁ * (3 - 32 * σ * (θ₂^2 * θ₃^6 + θ₂^6 * θ₃^2))))
            else
                error("Only order 0 and 1 are implemented for the multipole method.")
            end
        elseif network == "adjacent"
            if order == 0
                # Eq. 22 from Claesson and Javed 2019 (adjacent network, zeroth-order multipole)
                Ra = 2 * Rₚ + (1 / (π * kg)) * (log(s / ro) + σ * log((rb^2 + (s / 2)^2) /
                    (rb^2 - (s / 2)^2)))
            elseif order == 1
                # Eq. 23 from Claesson and Javed 2019 (adjacent network, first-order multipole)
                θ₁ = ro^2 / (4 * (s / 2)^2)
                θ₂ = (s / 2)^2 / (rb^8 - (s / 2)^8)^(1/4)
                θ₃ = rb^2 / (rb^8 - (s / 2)^8)^(1/4)
                b₁ = (1 - β) / (1 + β)
                V1 = 1 - 8 * σ * θ₂^3 * θ₃
                V2 = 3 + 8 * σ * θ₂ * θ₃^3
                M11 = 1 + 16 * b₁ * σ * θ₁ * (3 * θ₂^3 * θ₃^5 + θ₂^7 * θ₃)
                M21 = b₁ * θ₁
                M12 = -M21
                M22 = -1 - 16 * b₁ * σ * θ₁ * (θ₂ * θ₃^7 + 3 * θ₂^5 * θ₃^3)

                Ra = 2 * Rₚ + (1 / (π * kg)) * (log(s / ro) + σ * log((rb^2 + (s / 2)^2) /
                    (rb^2 - (s / 2)^2))) + (1 / (π * kg)) * ((b₁ * θ₁ * (V2^2 * M11 - 2 * V1 * V2 *
                    M21 - V1^2 * M22)) / (M11 * M22 + M21^2))
            else
                error("Only order 0 and 1 are implemented for the multipole method.")
            end
        else
            error("Only `diagonal` and `adjacent` networks exist in double U-loop configuration.")
        end
    else
        error("Only single and double U-loops are implemented for the multipole method.")
    end
    return Ra
end
function resistance_ULoop_total_internal(V::Real, s::Real, rb::Real, ro::Real, ri::Real,
    ks::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0;
    nLoop::Int=1, order::Int=1, network::String="diagonal")
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp)

    # Compute Ra with the first-order multipole method
    return resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop, order=order,
        network=network)
end

"""
    resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra)
    resistance_ULoop_effective(V, H, s, rb, ro, ks, kg, cf, ρf, Rp, Rf, nLoop=1)
    resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf, nLoop=1)

Function that computes the effective thermal borehole resistance (also named Rb*). Effective Rb
allows considering the thermal short-circuiting along the borehole. Two types of boundary 
conditions are commonly used: (1) uniform borehole wall temperature (UBW) or (2) uniform heat 
flux (UHF). The most practical approach is to use an average of both approach.
Note: The effective resistance formula is derived for single U-tube configurations (2 pipes per
borehole, `nLoop = 1`).
# Arguments
    - `V`: Fluid flow rate in pipe [m³/s]
    - `H`: Borehole length [m]
    - `s`: Shank spacing (distance between 2 legs of a U-tubes) [m]
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outlet radius [m]
    - `ri`: Pipe inner radius [m]
    - `ks`: Ground thermal conductivity [W/mK]
    - `kg`: Grout thermal conductivity [W/mK]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_μ(T)`)
    - `Rb`: Borehole thermal resistance [mK/W]
    - `Ra`: Total internal thermal resistance [mK/W]
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `nLoop`: Number of loops in the ground heat exchanger (2 pipes per loop) (default 1)
# Outputs
    - `Rbₑ`: Effective borehole thermal resistance [mK/W]
# Reference
    - Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal 
        resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
    - Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. 
        Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. 
        https://doi.org/10.1016/B978-0-08-100311-4.00003-0
"""
function resistance_ULoop_effective(V::Real, H::Real, cf::Real, ρf::Real, Rb::Real, Ra::Real)
    # UBW - See Eq. 3.68-3.70 of Javed et Spitler (2016)
    R1b = 2 * Rb                                    # Eq. 3.12
    R12 = (2 * Ra * R1b) / (2 * R1b - Ra)           # Eq. 3.14
    tmp = H / (V * cf * ρf)
    η = tmp / (2 * Rb) * sqrt(1 + 4 * Rb / R12)     # Eq. 3.69
    if η <= 1
        Rbₑ1 = Rb + (1 / (3 * R12)) * tmp^2 + (1 / (12 * Rb)) * tmp^2
    else
        Rbₑ1 = Rb * η * coth(η)
    end

    # UHF - See Eq. 3.67 of Javed et Spitler (2016)
    Rbₑ2 = Rb + (1 / (3 * Ra)) * tmp^2

    # Final calculation of the effective borehole thermal resistance
    return 0.5 * (Rbₑ1 + Rbₑ2)
end
function resistance_ULoop_effective(V::Real, H::Real, s::Real, rb::Real, ro::Real, ks::Real,
    kg::Real, cf::Real, ρf::Real, Rp::Real, Rf::Real; nLoop::Int=1)
    # Compute Rb and Ra with the first-order multipole methods
        Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop)
        Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop)

    # Compute Rbₑ
    return resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra)
end
function resistance_ULoop_effective(V::Real, H::Real, s::Real, rb::Real, ro::Real, ri::Real,
    ks::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0;
    nLoop::Int=1)
    # Compute fluid and pipe resistances
    Rf = resistance_fluid(V / (π * ri^2), ri, kf, cf, ρf, μf, ϵ)
    Rp = resistance_pipe(ro, ri, kp)

    # Compute Rb and Ra with the first-order multipole method
    Rb = resistance_ULoop_borehole(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop)
    Ra = resistance_ULoop_total_internal(s, rb, ro, ks, kg, Rp, Rf; nLoop=nLoop)

    # Compute Rbₑ
    return resistance_ULoop_effective(V, H, cf, ρf, Rb, Ra)
end

"""
    resistance_coaxial(rii, rio, roi, roo, rb, kg, kpi, kpo, hin, hann)
    resistance_coaxial(V, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ=0.0)

Computes the two thermal resistances of the resistance network of a coaxial (concentric-tube)
ground heat exchanger, following the model of Lamarche (2021):
- `R12`: internal resistance between the center pipe and the annulus (Eq. 1);
- `R1` : resistance between the annulus fluid and the borehole wall (Eq. 2).
The center pipe carries the fluid inside the inner pipe, while the annulus is the region between
the inner pipe outer wall and the outer pipe inner wall. Following Eq. 8 of Lamarche (2021), `R1`
is also the (steady) borehole resistance `Rb` of a coaxial exchanger. Both flow directions
("center-in" and "annulus-in") share the same `R1`/`R12` network.
# Arguments
    - `V`: Fluid flow rate in the exchanger [m³/s]
    - `rii`: Inner pipe inner radius [m]
    - `rio`: Inner pipe outer radius [m]
    - `roi`: Outer pipe inner radius [m]
    - `roo`: Outer pipe outer radius [m]
    - `rb`: Borehole radius (`rb > roo`, includes grout thickness) [m]
    - `kg`: Grout thermal conductivity (if present) [W/mK]
    - `kpi`: Inner pipe thermal conductivity [W/mK]
    - `kpo`: Outer pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `hin`: Convection coefficient inside the inner (center) pipe [W/m²K]
    - `hann`: Convection coefficient in the annulus region [W/m²K]
# Output
    - `R1`: Annulus-to-borehole-wall resistance (= coaxial borehole resistance `Rb`) [mK/W]
    - `R12`: Center-to-annulus internal resistance [mK/W]
# Reference
    - Lamarche, L. (2021). Analytic models and effective resistances for coaxial ground heat
        exchangers. Geothermics, 97, 102224. https://doi.org/10.1016/j.geothermics.2021.102224
"""
function resistance_coaxial(rii::Real, rio::Real, roi::Real, roo::Real, rb::Real, kg::Real,
    kpi::Real, kpo::Real, hin::Real, hann::Real)
    # Center-to-annulus internal resistance (Eq. 1 of Lamarche 2021):
    # convection inside inner pipe + inner pipe wall conduction + annulus convection (inner face)
    R12 = 1 / (hin * 2 * π * rii) + resistance_pipe(rio, rii, kpi) + 1 / (hann * 2 * π * rio)

    # Annulus-to-borehole-wall resistance (Eq. 2 of Lamarche 2021):
    # annulus convection (outer face) + outer pipe wall conduction + grout conduction
    R1 = 1 / (hann * 2 * π * roi) + resistance_pipe(roo, roi, kpo) + resistance_pipe(rb, roo, kg)

    return R1, R12
end
function resistance_coaxial(V::Real, rii::Real, rio::Real, roi::Real, roo::Real, rb::Real,
    kg::Real, kpi::Real, kpo::Real, kf::Real, cf::Real, ρf::Real, μf::Real, ϵ::Real=0.0)
    # Center pipe convection coefficient (circular pipe, hydraulic diameter = 2·rii)
    V̇in = V / (π * rii^2)                        # Mean fluid speed in the center pipe [m/s]
    Nu_in = Nusselt(V̇in, rii, kf, cf, ρf, μf, ϵ)
    hin = convection_coefficient(Nu_in, rii, kf)

    # Annulus convection coefficient (hydraulic diameter Dₕ = 2·(roi - rio))
    V̇ann = V / (π * (roi^2 - rio^2))             # Mean fluid speed in the annulus [m/s]
    Nu_ann = Nusselt_annulus(V̇ann, roi, rio, kf, cf, ρf, μf, ϵ, ϵ)
    hann = convection_coefficient(Nu_ann, roi - rio, kf)  # Same h applied to both annulus faces

    return resistance_coaxial(rii, rio, roi, roo, rb, kg, kpi, kpo, hin, hann)
end

"""
    resistance_coaxial_effective(V, H, cf, ρf, R1, R12; model="UHF")
    resistance_coaxial_effective(V, H, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf,
        μf, ϵ=0.0; model="UHF")

Computes the effective borehole thermal resistance (`Rb*`) of a coaxial ground heat exchanger,
accounting for the axial thermal short-circuit between the center pipe and the annulus, following
Lamarche (2021). Four closed-form models are available through the `model` keyword:
- `"UHF"` (uniform heat flux, Eq. 31) — the recommended compromise for coaxial exchangers, valid
    for a uniform far-field temperature;
- `"UBW"` (uniform borehole wall temperature, Eq. 14) - not as accurate as `"UHF"`;
- `"mean"` — the average of both, matching the convention of `resistance_ULoop_effective`;
- `"UHF_gradient"` (linearly-varying heat flux, Eqs. 58-63) — refines `"UHF"` when the far-field
    temperature itself varies linearly with depth (e.g. a geothermal gradient).
Both flow configurations ("center-in" and "annulus-in") yield the same effective resistance under
the `"UHF"`/`"UBW"`/`"mean"` steady-flux assumptions (Section 2.2 of Lamarche 2021); `"UHF_gradient"`
is the exception, by construction (its whole point is to distinguish the two flow directions under
a gradient).
# Arguments
    - `V`: Fluid flow rate in the exchanger [m³/s]
    - `H`: Borehole length [m]
    - `rii`: Inner pipe inner radius [m]
    - `rio`: Inner pipe outer radius [m]
    - `roi`: Outer pipe inner radius [m]
    - `roo`: Outer pipe outer radius [m]
    - `rb`: Borehole radius (`rb > roo`, includes grout thickness) [m]
    - `kg`: Grout thermal conductivity (if present) [W/mK]
    - `kpi`: Inner pipe thermal conductivity [W/mK]
    - `kpo`: Outer pipe thermal conductivity [W/mK]
    - `kf`: Fluid thermal conductivity [W/mK] (`water_k(T)`)
    - `cf`: Fluid specific heat [J/kgK] (`water_cp(T)`)
    - `ρf`: Fluid density [kg/m³] (`water_ρ(T)`)
    - `μf`: Fluid dynamic viscosity [kg/m/s] (`water_μ(T)`)
    - `ϵ`: Pipe roughness [m] (default 0.0)
    - `R1`: Annulus-to-borehole-wall resistance (coaxial `Rb`) [mK/W]
    - `R12`: Center-to-annulus internal resistance [mK/W]
    - `model`: Effective resistance model, `"UHF"`, `"UBW"`, `"mean"` or `"UHF_gradient"`
        (default `"UHF"`)
# Output
    - `Rbₑ`: Effective borehole thermal resistance [mK/W]
# Reference
    - Lamarche, L. (2021). Analytic models and effective resistances for coaxial ground heat
        exchangers. Geothermics, 97, 102224. https://doi.org/10.1016/j.geothermics.2021.102224
"""
function resistance_coaxial_effective(V::Real, H::Real, cf::Real, ρf::Real, R1::Real,
    R12::Real; model::String="UHF")
    # Dimensionless groups (Eq. 7 of Lamarche 2021)
    γ = H / (2 * ṁcf * R1)
    Ra = 4 * R1 * R12 / (4 * R1 + R12)
    ξ = sqrt(Ra / (4 * R1))
    η = γ / ξ

    if model == "UBW"                   # Uniform borehole wall temperature (Eq. 14)
        return R1 * η * coth(η)
    elseif model == "UHF"               # Uniform heat flux (Eq. 31)
        return R1 * (1 + (Ra / R12) * η^2 / 3)
    elseif model == "mean"              # Average of both models
        return 0.5 * (R1 * η * coth(η) + R1 * (1 + (Ra / R12) * η^2 / 3))
    elseif model == "UHF_gradient"      # Linear heat flux under a geothermal gradient (Eqs. 58-63)
        return R1 * (1 + H / (6 * V * ρf * cf   * R1) + H^2 / (4 * (V * ρf * cf  )^2 * R1 * R12))
    else
        error("Only `UHF`, `UBW`, `mean` and `UHF_gradient` models are implemented for coaxial " *
            "exchangers.")
    end
end
function resistance_coaxial_effective(V::Real, H::Real, rii::Real, rio::Real, roi::Real,
    roo::Real, rb::Real, kg::Real, kpi::Real, kpo::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    ϵ::Real=0.0; model::String="UHF")
    # Compute the coaxial R1 and R12 resistances
    R1, R12 = resistance_coaxial(V, rii, rio, roi, roo, rb, kg, kpi, kpo, kf, cf, ρf, μf, ϵ)

    # Compute Rbₑ
    return resistance_coaxial_effective(V, H, cf, ρf, R1, R12; model=model)
end