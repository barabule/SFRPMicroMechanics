module SFRPMicroMechanics

using LinearAlgebra
using Statistics
using GeometryBasics
using StaticArrays


# Include sub-modules
include("typedefs.jl")
include("eshelby.jl")
include("thermal_expansion.jl")
include("utils.jl")
include("closures.jl")
include("homogenization.jl")
include("analyse.jl")


# Re-export main functions for the user
export compute_orthotropic_properties, 
       calibrate_mat215, 
       simulate_tensile_test,
       generate_stiffness_mesh

end # module SFRPMicroMechanics
