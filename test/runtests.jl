 import Pkg; Pkg.activate("test")
 using Test
 using SFRPMicroMechanics

@testset "ALL TESTS" verbose = true begin
include("testsets.jl")
end