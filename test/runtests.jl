 import Pkg; Pkg.activate("test")
 using Test
 using SFRPMicroMechanics
 using LinearAlgebra, Random

@testset "ALL TESTS" verbose = true begin
include("testsets.jl")
end