using CoolProp

const P = 100000 # [pa]

const _MIXTURE_PREFIX = Dict{Symbol,String}(
    :MPG => "MPG", # Propylene glycol / water
    :MEG => "MEG", # Ethylene glycol / water
    :MMA => "MMA", # Methanol / water
    :MEA => "MEA", # Ethanol / water
    :MKA => "MKA", # Potassium acetate / water
    :MKF => "MKF", # Potassium formate / water
)

"""
    fluid_property(T, fluid; percentage)

Function that returns the four properties of the heat transfer fluid (thermal conductivity, specific
heat capacity, density and viscosity), given the temperature, the type of the fluid (pure or mix),
and the percentage (if the fluid is a mixture). The fluid can be pure water (:water), or a mix
of propylene glycol (:MPG), ethylene glycol (:MEG), methanol (:MMA), ethanol (:MEA), potassium
acetate (:MKA), or potassium formate (:MKF).
# Arguments
    - `T`: Fluid's temperature [°C]
    - `fluid`: Fluid's symbol (pure: :water) or (mix: :MPG, :MEG, :MMA, :MEA, :MKA, :MKF) [-]
    - `percentage`: Fluid's concentration, mass fraction of the additive [%m], default of 100%
      (ignored for :water)
# Output
    - `k`: Fluid thermal conductivity [W/mK]
    - `cp`: Fluid specific heat capacity [J/kgK]
    - `ρ`: Fluid's density [kg/m³]
    - `μ`: Fluid's viscosity [kg/ms]
# Reference
    - CoolProp, Incompressible fluids, 2026.
        https://coolprop.org/fluid_properties/Incompressibles.html
# Example
    ```julia
    julia> fluid_property(20, :water)
    (0.598, 4181.3, 998.2, 0.001002)

    julia> fluid_property(20, :MPG; percentage=50)
    (0.406, 3600.0, 1050.0, 0.0045)
    ```
"""
function fluid_property(T::Real, fluid::Symbol=:water; percentage::Real = 100)
    # Check that the temperature is within the valid range for CoolProp
    T_K = T + 273.15 # CoolProp expects kelvin, but this package works in °C throughout.

    # Check that the percentage is within the valid range for CoolProp
    fluid_string = if fluid == :water
        T_boil = PropsSI("T", "P", P, "Q", 0, "Water") - 273.15 # [°C] boiling point (99.61°C 1 bar)
        if T < 0 || T >= T_boil
            @warn "Temperature out of the liquid range (0°C ≤ T < $(round(T_boil, digits=1))°C "
        end
        "Water"
    elseif haskey(_MIXTURE_PREFIX, fluid)
        "INCOMP::" * _MIXTURE_PREFIX[fluid] * "-" * string(percentage) * "%"
    else
        error("Unsupported fluid :$fluid. Supported fluids are :water, " *
            join(sort(":" .* string.(keys(_MIXTURE_PREFIX))), ", ") * ".")
    end

    k = PropsSI("L", "P", P, "T", T_K, fluid_string)
    cp = PropsSI("C", "P", P, "T", T_K, fluid_string)
    ρ = PropsSI("D", "P", P, "T", T_K, fluid_string)
    μ = PropsSI("V", "P", P, "T", T_K, fluid_string)

    return k, cp, ρ, μ
end

"""
    water_k(T::Real)

!!! warning "Deprecated"
    Use [`fluid_property`](@ref)`(T, :water)` instead.

Water thermal conductivity k(T), 0 ≤ T ≤ 99.6°C at 1 bar. The polynomial equation is a fit to the
data from the Engineering Toolbox. The values can be validates with a temperature vector
`T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - k: Thermal conductivity [W/mK]
# Reference
    - The Engineering ToolBox (2018). Thermal Conductivity of Water: Temperature and Pressure Data.
        [online] Available at: https://www.engineeringtoolbox.com/water-liquid-gas-thermal-
        conductivity-temperature-pressure-d_2012.html [Accessed 2026-01-14].
"""
function water_k(T::Real)
    Base.depwarn("`water_k` is deprecated, use `fluid_property(T, :water)` instead.", :water_k)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_k` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 0.5557250521318174 + 0.002490814640452007*T - 2.1170044416971473e-5*T^2 +
        1.285515973680875e-7*T^3 - 4.546428806458628e-10*T^4 + 9.750314739837196e-14*T^5
end

"""
    water_cp(T::Real)

!!! warning "Deprecated"
    Use [`fluid_property`](@ref)`(T, :water)` instead.

Water specific heat capacity cp(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the isobaric
specific heat capacity data from the Engineering Toolbox. The values can be validates with a
temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - cp: Specific heat capacity [J/(kg·K)]
# Reference
    - The Engineering ToolBox (2004). Specific Heat Capacity of Water: Temperature-Dependent Data
        and Calculator. [online] Available at:
        https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
        [Accessed 2026-01-14].
"""
function water_cp(T::Real)
    Base.depwarn("`water_cp` is deprecated, use `fluid_property(T, :water)` instead.", :water_cp)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_cp` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 4219.849078078278 - 3.266686616602623*T + 0.09969277880041719*T^2 -
        0.0014860911377001344*T^3 + 1.161963811563561e-5*T^4 - 3.5034316470844105e-8*T^5
end

"""
    water_ρ(T::Real)

!!! warning "Deprecated"
    Use [`fluid_property`](@ref)`(T, :water)` instead.

Water density ρ(T), 0 ≤ T ≤ 100°C at 1 atm. The polynomial equation is a fit to the data from the
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - ρ: Density [kg/m³]
# Reference
    - The Engineering ToolBox (2003). Water Density, Specific Weight and Thermal Expansion
        Coefficients - Temperature and Pressure Dependence. [online] Available at:
        https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
        [Accessed 2026-01-14].
"""
function water_ρ(T::Real)::Float64
    Base.depwarn("`water_ρ` is deprecated, use `fluid_property(T, :water)` instead.", :water_ρ)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_ρ` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 999.8475436930158 + 0.06180756931966996*T - 0.008309049138917115*T^2 +
        6.35713412865478e-5*T^3 - 3.8404497053894326e-7*T^4 + 1.0249871031879443e-9*T^5
end

"""
    water_μ(T::Real)

!!! warning "Deprecated"
    Use [`fluid_property`](@ref)`(T, :water)` instead.

Dynamic viscosity μ(T), 0 ≤ T ≤ 100°C. The polynomial equation is a fit to the data from the
Engineering Toolbox. The values can be validates with a temperature vector `T = 0.1:1:100`.
# Argument
    - T: Temperature [°C]
# Output
    - μ: Dynamic viscosity [Pa·s or kg/(m·s)]
# Reference
    - The Engineering ToolBox (2004). Water - Dynamic and Kinematic Viscosity at Various
        Temperatures and Pressures. [online] Available at: https://www.engineeringtoolbox.com/water-
        dynamic-kinematic-viscosity-d_596.html [Accessed 2026-01-14].
"""
function water_μ(T::Real)
    Base.depwarn("`water_μ` is deprecated, use `fluid_property(T, :water)` instead.", :water_μ)
    if T < 0 || T > 100
         @warn "Temperature out of range (0 ≤ T ≤ 100°C).
         `water_μ` fits data from Engineering Toolbox, and may not be accurate outside this range."
    end
    return 0.001790966556989398 - 5.965082369793418e-5*T + 1.3185191782991122e-6*T^2 -
        1.8236868027209892e-8*T^3 + 1.3644271817518522e-10*T^4 - 4.137645533574321e-13*T^5
end
