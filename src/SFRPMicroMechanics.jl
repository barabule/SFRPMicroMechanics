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
export IsotropicElasticParameters, OrthotropicElasticParameters, stiffness_matrix_voigt, extract_orthotropic_constants
export is_isotropic, compute_orthotropic_properties
export eshelby_tensor, SphericalInclusion, SpheroidalInclusion, DiscInclusion, ThinDiscInclusion, NeedleInclusion
export ThermalExpansion, halpin_tsai, orientation_average, mori_tanaka

end # module SFRPMicroMechanics
