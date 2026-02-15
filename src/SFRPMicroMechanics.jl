module SFRPMicroMechanics

using LinearAlgebra
using Statistics
using GeometryBasics
using StaticArrays
using Tullio

# Include sub-modules
include("homogenization.jl")


# Re-export main functions for the user
export compute_orthotropic_properties, 
       calibrate_mat215, 
       simulate_tensile_test,
       generate_stiffness_mesh

end # module SFRPMicroMechanics
