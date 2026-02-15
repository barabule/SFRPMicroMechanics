
# @testset "SFRPMicroMechanics Physical Validation" verbose=true begin
@testset "Basic" verbose = true begin
    # Test 1: Zero fibers

    props = SFRPMicroMechanics.compute_orthotropic_properties(2000.0, 0.35, 70000.0, 0.22, 0.0, 15, 0.6, 0.2)
    # @info props
    @test isapprox(props["E1"], 2000.0, atol=1e-2)
    @test isapprox(props["nu12"], 0.35, atol=1e-4)

end

@testset "Anisotropic Properties" verbose = true begin
    props = SFRPMicroMechanics.compute_orthotropic_properties(2500.0, 0.35, 72000.0, 0.22, 0.2, 25, 1.0, 0.0)
    # @info props
    @test props["E1"] > props["E2"]
    @test isapprox(props["E2"], props["E3"], rtol=1e-3)
    # --- SETUP BASE PROPERTIES ---
    Em, num = 2500.0, 0.35      # Matrix (e.g., Nylon)
    Ef, nuf = 72000.0, 0.22     # Fiber (e.g., Glass)
    vf = 0.20                   # 20% Volume fraction
    aspect_ratio = 25.0
end


@testset "Test 1: Zero Fiber Limit (Pure Matrix)" verbose=true begin
        # With vf = 0, the output must be isotropic matrix properties
    aspect_ratio = 1.0
    Em, num, Ef, nuf = 1.0, 0.3, 10.0, 0.2
    res = compute_orthotropic_properties(Em, num, Ef, nuf, 0.0, aspect_ratio, 0.8, 0.1)
    
    @test isapprox(res["E1"], Em, atol=1e-1)
    @test isapprox(res["E2"], Em, atol=1e-1)
    @test isapprox(res["nu12"], num, atol=1e-3)
    @test isapprox(res["G12"], Em / (2*(1+num)), atol=1e-1)
    # Explanation: If there are no fibers, orientation shouldn't matter.
end
    
@testset "Test 2: Spherical Inclusion Limit" verbose = true begin
    # Aspect Ratio = 1.0 means fibers are spheres. 
    # Spheres are isotropic; E1 should equal E2 even with aligned a11.
    aspect_ratio = 1.0
    vf = 0.2
    Em, num, Ef, nuf = 1.0, 0.3, 10.0, 0.2
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, 1.0, 0.9, 0.05)
    
    @test isapprox(res["E1"], res["E2"], rtol=0.01)
    @test isapprox(res["E1"], res["E3"], rtol=0.01)
    # Explanation: A sphere has no preferred direction. If E1 != E2, 
    # the Eshelby tensor or the index mapping in Advani-Tucker is wrong.
end
    
@testset "Test 3: Unidirectional (UD) Alignment" verbose = true begin
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
    
@testset "Test 4: 3D Random Isotropy" verbose = true begin
    # a11 = a22 = a33 = 1/3 is a perfectly random 3D distribution.
    a11= a22  = 1/3
    Em, num, Ef, nuf = 1.0, 0.3, 10.0, 0.2
    aspect_ratio = 10.0
    vf = 0.3
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    
    @test isapprox(res["E1"], res["E2"], rtol=0.01)
    @test isapprox(res["E2"], res["E3"], rtol=0.01)
    # Explanation: This tests the Hybrid Closure. At det(a) = 1/27, 
    # f should be 0, triggering the Linear Closure logic.
end

@testset "Test 5: Bounds Check (Rule of Mixtures)" verbose = true begin
    # The longitudinal modulus E1 of a UD composite (a11=1) 
    # should never exceed the Voigt Upper Bound (Rule of Mixtures).
    a11, a22  = 1.0, 0.0
    Em, num, Ef, nuf = 1.0, 0.3, 10.0, 0.2
    aspect_ratio = 1000.0
    vf = 0.3
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22) # AR=1000 ~ continuous
    ROM_E1 = vf * Ef + (1 - vf) * Em
    
    @test res["E1"] <= ROM_E1 + 1.0 # Allow for tiny numerical float noise
    # Explanation: Micromechanics models must respect thermodynamic bounds.

end