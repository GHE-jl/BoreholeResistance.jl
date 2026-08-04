using CoolProp

const P = 100000 # [pa]
"""
    fluid_property(T, fluid; percentage)

Function that returns the four properties of the fluid (thermal conductivity, volumetric heat 
capacity, density and viscosity), given the temperature, the type of the fluid (pure or mix), and 
the percentage (if the fluid is a mixture). The fluid can be pure water (:water), or an aqueous mix 
of propylene glycol (:MPG), or ethylene glycol (:MEG) .
# Arguments
    - `T`: Fluid's temperature [K]
    - `fluid`: Fluid's symbole (pure: :water) or (mix: :MPG, :MEG) [-]
    - `percentage`: Fluid's concentration [%m], default of 100% if not specified
# Output
    - `k`: Fluid thermal conductivity [W/mK]
    - `cs`: Fluid volumetric heat capacity [J/Km³]
    - `ρ`: Fluid's density [kg/m³]
    - `μ`: Fluid's viscosity [kg/ms] 
# Validity ranges (CoolProp ranges)
    - `T`:
        - :water : 273.15     - 473.15 [K]
        - :MPG :   252.581755 - 373.15 [K]
        - :MEG :   223.147354 - 373.15 [K]
    -`Percentage`:
        - :water : 100 [%m]
        - :MPG : 0-60 [%m]
        - :MEG : 0-60 [%m]
    
# Reference
    - CoolProp, Incompresible fluids, 2026.
    https://coolprop.org/fluid_properties/Incompressibles.html
"""
function fluid_property(T::Real, fluid::Symbol; percentage::Int = 100)
    if fluid == :water
        (273.15 <= T <= 473.15) ? 
        fluid = "Water" : 
        @warn("Water's temperature is outside the permitted range of 273.15-473.15K")

    elseif fluid == :MPG
        ((252.581755 <= T <= 373.15) & (0 <= percentage <= 60)) ? 
        fluid = "INCOMP::MPG-"*string(percentage)*"%" : 
        @warn("MPG's input are outside the permitted range of 252.581755-373.15K and 0-60%") & 
        return
    
    elseif fluid == :MEG
        ((223.147354 <= T <= 373.15) & (0 <= percentage <= 60)) ?
        fluid = "INCOMP::MEG-"*string(percentage)*"%" : 
        @warn("MEG's input are outside the permitted range of 223.147354-373.15K and 0-60%") &
        return
    
    else
        @warn("This fluid is not yet supported, water will be used instead")
        fluid = "Water"

    end

    k = PropsSI("L", "P", P, "T", T, fluid)
    cs = PropsSI("C", "P", P, "T", T, fluid)
    ρ = PropsSI("D", "P", P, "T", T, fluid)
    μ = PropsSI("V", "P", P, "T", T, fluid)
    
    return k, cs, ρ, μ
end