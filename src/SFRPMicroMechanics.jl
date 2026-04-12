module SFRPMicroMechanics

using LinearAlgebra
using Statistics
using GeometryBasics
using StaticArrays
using Random
using Tensors
using Combinatorics: permutations


# Include sub-modules
include("types_stiffness.jl")
include("types_thermal_expansion.jl")
include("orientation_tensors.jl")
include("eshelby.jl")

include("constituent.jl")

include("utils.jl")
include("closures.jl")
include("homogenization.jl")
include("thermal_expansion.jl")

include("analyse.jl")


# Re-export main functions for the user
export mori_tanaka
export FiberPhase
#elastic 
export IsotropicProperties, OrthotropicProperties, TransverseIsotropicProperties
export stiffness_matrix_voigt, extract_orthotropic_constants

#inclusions
export eshelby_tensor, SphericalInclusion, SpheroidalInclusion
export DiscInclusion, ThinDiscInclusion, NeedleInclusion

export ThermalExpansion, compute_all_thermal_expansions

export halpin_tsai, orientation_average
#util
export apparent_modulus

#closures
export LinearClosure, QuadraticClosure, HybridClosure, HL1Closure, HL2Closure
export ORF, ORS, ORL, ORFM, ORW, ORW3, IBOF


export PrincipalOrientationTensor, FullOrientationTensor, decompose


end # module SFRPMicroMechanics
