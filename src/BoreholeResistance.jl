module BoreholeResistance

# Fluid convective resistance (Reynolds, Prandtl, Nusselt, friction factor)
include("resistance_fluid.jl")
# Pipe conduction resistance
include("resistance_pipe.jl")
# Borehole thermal resistance (multipole method) and effective resistance
include("resistance_borehole.jl")

# Fluid properties and pipe flow utilities
include("utils.jl")

# Fluid resistance exports
export Reynolds, Prandtl, Nusselt, Nusselt_annulus,
    friction_factor_Colebrook_White, friction_factor_Tkachenko_Mileikovskyi,
    resistance_fluid

# Pipe resistance exports
export resistance_pipe

# Borehole resistance exports
export resistance_borehole_multipole, resistance_total_internal_multipole,
    resistance_borehole_effective

# Utility exports
export water_ρ, water_cp, water_k, water_μ, head_loss_Darcy_Weisbach

end