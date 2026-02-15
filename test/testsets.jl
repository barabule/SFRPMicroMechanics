
# @testset "SFRPMicroMechanics Physical Validation" verbose=true begin
@testset "No Fibers" verbose = true begin
    # Test 1: Zero fibers
    #vf = 0
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.0, 15.0, 0.7, 0.25
    props = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    # @info props
    @test isapprox(props["E1"], Em, atol=1e-2)
    @test isapprox(props["nu12"], num, atol=1e-4)
    @test isapprox(props["E2"], Em, atol=1e-1)
    @test isapprox(props["G12"], Em / (2*(1+num)), atol=1e-1)
end

@testset "Anisotropic Properties" verbose = true begin
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.3, 15.0, 0.7, 0.25
    props = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    # @info props
    @test props["E1"] > props["E2"]
    @test isapprox(props["E2"], props["E3"], rtol=1e-3)
    
    
end

    
@testset "Spherical Inclusion Limit" verbose = true begin
    # Aspect Ratio = 1.0 means fibers are spheres. 
    # Spheres are isotropic; E1 should equal E2 even with aligned a11.
    ar = 1.0
    vf = 0.2
    Em, num, Ef, nuf = 2000, 0.3, 70e3, 0.22
    a11 = 0.7
    a22 = 0.25
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    for (k,v) in res
        @info "$k = $v"
    end
    @test isapprox(res["E1"], res["E2"], rtol=0.01)
    @test isapprox(res["E1"], res["E3"], rtol=0.01)
    # Explanation: A sphere has no preferred direction. If E1 != E2, 
    # the Eshelby tensor or the index mapping in Advani-Tucker is wrong.
end
    
@testset "Unidirectional (UD) Alignment" verbose = true begin
    # a11 = 1.0 means all fibers are perfectly parallel.
    # This must yield a Transversely Isotropic material.
    a11, a22  = 1.0,  0.0
    Em, num, Ef, nuf = 1.0, 0.3, 10.0, 0.2
    aspect_ratio = 10.0
    vf = 0.3
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    
    @test res["E1"] > res["E2"]
    @test isapprox(res["E2"], res["E3"], rtol=1e-4)
    @test isapprox(res["nu12"], res["nu13"], rtol=1e-4)
    # Explanation: In a UD composite, the properties perpendicular 
    # to the fiber (2 and 3 directions) must be identical.
end
    
@testset "3D Random Isotropy" verbose = true begin
    # a11 = a22 = a33 = 1/3 is a perfectly random 3D distribution.
    a11= a22  = 1/3
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 10.0
    vf = 0.3
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    
    @test isapprox(res["E1"], res["E2"], rtol=0.01)
    @test isapprox(res["E2"], res["E3"], rtol=0.01)
    # Explanation: This tests the Hybrid Closure. At det(a) = 1/27, 
    # f should be 0, triggering the Linear Closure logic.
end

@testset "Bounds Check (Rule of Mixtures)" verbose = true begin
    # The longitudinal modulus E1 of a UD composite (a11=1) 
    # should never exceed the Voigt Upper Bound (Rule of Mixtures).
    a11, a22  = 1.0, 0.0
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 1000.0
    vf = 0.3
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22) # AR=1000 ~ continuous
    ROM_E1 = vf * Ef + (1 - vf) * Em
    @info "E1 computed = $(res["E1"])"
    @info "Upper bound = $ROM_E1"
    @test res["E1"] <= ROM_E1 + 1.0 # Allow for tiny numerical float noise
    # Explanation: Micromechanics models must respect thermodynamic bounds.

end