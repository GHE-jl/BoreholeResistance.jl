module BoreholeResistance

# Fluid convective resistance (Reynolds, Prandtl, Nusselt, friction factor)
include("resistance_fluid.jl")
# Pipe conduction resistance
include("resistance_pipe.jl")
# Borehole thermal resistance (multipole method) and effective resistance
include("resistance_borehole.jl")

# Fluid properties from CoolProp (water and antifreeze mixtures), plus legacy deprecated
# polynomial fits (water_k, water_cp, water_ρ, water_μ)
include("fluid_property.jl")

# Fluid resistance exports
export Reynolds, Prandtl, Nusselt, Nusselt_annulus,
    friction_factor_Colebrook_White, friction_factor_Tkachenko_Mileikovskyi, convection_coefficient,
    resistance_fluid

# Pipe resistance exports
export resistance_pipe

# Borehole resistance exports
export resistance_ULoop_borehole, resistance_ULoop_total_internal,
    resistance_ULoop_effective,
    resistance_coaxial, resistance_coaxial_effective

# Fluid property exports
export fluid_property
# Deprecated, use fluid_property instead
export water_ρ, water_cp, water_k, water_μ

end