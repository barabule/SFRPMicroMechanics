@testset "Tensor Conversions" verbose = true begin
    C = SFRPMicroMechanics.isotropic_stiffness(210.0, 0.3)
    tens = SFRPMicroMechanics.convert_66_to_3333(C; mandel = false)
    Cc = SFRPMicroMechanics.convert_3333_to_66(tens;mandel = false)
    tensc = SFRPMicroMechanics.convert_66_to_3333(Cc; mandel = false)

    @test all(C .≈ Cc)
    @test all(tens .≈ tensc)


end

# @testset "SFRPMicroMechanics Physical Validation" verbose=true begin
@testset "No Fibers" verbose = true begin
    # Test 1: Zero fibers
    #vf = 0
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.0, 15.0, 0.7, 0.25
    props = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    # @info props
    @test isapprox(props.E1, Em, atol=1e-2)
    @test isapprox(props.nu21, num, atol=1e-4)
    @test isapprox(props.E2, Em, atol=1e-1)
    @test isapprox(props.G12, Em / (2*(1+num)), atol=1e-1)
end

@testset "Eshelby Tensor" verbose = true begin
    AR = 1e6 #S1111 -> 0 when AR-> Inf
    S  = SFRPMicroMechanics.eshelby_tensor_spheroid(0.3, AR)
    @test S[1,1] <= sqrt(eps(Float64))
    
    AR = 1e-10 #S1111 -> 1 when AR-> 0
    S = SFRPMicroMechanics.eshelby_tensor_spheroid(0.3, AR)
    @test S[1,1] ≈ 1
end


@testset "Anisotropic Properties" verbose = true begin
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.3, 15.0, 0.7, 0.25
    props = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    # @info props
    @test props.E1 > props.E2
    @test isapprox(props.E2, props.E3, rtol=1e-3)
    
    
end

@testset "Extracted Props" verbose = true begin
    p = SFRPMicroMechanics.OrthotropicElasticParameters(;E1 = 10.0,
                                                        E2 = 1.0,
                                                        E3 = 0.9,
                                                        G12 = 0.5,
                                                        G23 = 0.7,
                                                        G31 = 0.3,
                                                        nu21 = 0.05,
                                                        nu32 = 0.2,
                                                        nu31 = 0.1)
    Cm = SFRPMicroMechanics.stiffness_matrix_voigt(p)
    pex = SFRPMicroMechanics.extract_orthotropic_constants(Cm)

    @test p.E1 ≈ pex.E1
    @test p.E2 ≈ pex.E2
    @test p.E3 ≈ pex.E3
    @test p.G12 ≈ pex.G12
    @test p.G23 ≈ pex.G23
    @test p.G31 ≈ pex.G31
    @test p.nu21 ≈ pex.nu21
    @test p.nu31 ≈ pex.nu31
    @test p.nu32 ≈ pex.nu32

end


    
@testset "Spherical Inclusion Limit" verbose = true begin
    # Aspect Ratio = 1.0 means fibers are spheres. 
    # Spheres are isotropic; E1 should equal E2 even with aligned a11.
    ar = 1.0
    vf = 0.15
    Em, num, Ef, nuf = 2000, 0.3, 70e3, 0.22
    a11 = 0.7
    a22 = 0.25
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    @info res
    @test isapprox(res.E1, res.E2, rtol=0.01)
    @test isapprox(res.E1, res.E3, rtol=0.01)
    # Explanation: A sphere has no preferred direction. If E1 != E2, 
    # the Eshelby tensor or the index mapping in Advani-Tucker is wrong.
end
    
@testset "Unidirectional (UD) Alignment" verbose = true begin
    # a11 = 1.0 means all fibers are perfectly parallel.
    # This must yield a Transversely Isotropic material.
    a11, a22  = 1.0,  0.0
    Em, num, Ef, nuf = 2500.0, 0.3, 70e3, 0.2
    aspect_ratio = 10.0
    vf = 0.17
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    
    @test res.E1 > res.E2
    @test isapprox(res.E2, res.E3, rtol=1e-4)
    @test isapprox(res.nu21, res.nu31, rtol=1e-4)
    # Explanation: In a UD composite, the properties perpendicular 
    # to the fiber (2 and 3 directions) must be identical.
end
    
@testset "3D Random Isotropy" verbose = true begin
    # a11 = a22 = a33 = 1/3 is a perfectly random 3D distribution.
    a11 = a22  = 1/3
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 10.0
    vf = 0.2
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    
    @test isapprox(res.E1, res.E2, rtol=0.01)
    @test isapprox(res.E2, res.E3, rtol=0.01)
    # Explanation: This tests the Hybrid Closure. At det(a) = 1/27, 
    # f should be 0, triggering the Linear Closure logic.
end

@testset "Bounds Check (Rule of Mixtures)" verbose = true begin
    # The longitudinal modulus E1 of a UD composite (a11=1) 
    # should never exceed the Voigt Upper Bound (Rule of Mixtures).
    a11, a22  = 1.0, 0.0
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 1000.0
    vf = 0.2
    res = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22) # AR=1000 ~ continuous
    ROM_E1 = vf * Ef + (1 - vf) * Em
    @info "E1 computed = $(res.E1)"
    @info "Upper bound = $ROM_E1"
    @test res.E1 <= ROM_E1 + 1.0 # Allow for tiny numerical float noise
    # Explanation: Micromechanics models must respect thermodynamic bounds.

end

@testset "Apparent Elastic Modulus" verbose = true begin

    p = SFRPMicroMechanics.OrthotropicElasticParameters(;E1 = 10.0,
                                                        E2 = 1.0,
                                                        E3 = 0.9,
                                                        G12 = 0.5,
                                                        G23 = 0.7,
                                                        G31 = 0.3,
                                                        nu21 = 0.05,
                                                        nu32 = 0.2,
                                                        nu31 = 0.1)

    angle = 0.0
    @test SFRPMicroMechanics.apparent_modulus(angle, p) ≈ p.E1
    angle = 90.0
    @test SFRPMicroMechanics.apparent_modulus(angle, p) ≈ p.E2

    #from LS-Dyna run with mat_002 with the same properties
    #matches somewhat for larger angles (close to 90)
    ref = Dict(0 => 2.82 / 20 / 0.01453,
                10=> 2.04 / 20 / 0.01453,
                20 => 1.09 / 20 / 0.01453,
                30 => 0.64 / 20 / 0.01453,
                40 => 0.445 / 20 / 0.01453,
                50 => 0.355 / 20 / 0.01453,
                60 => 0.312 / 20 / 0.01453,
                70 => 0.294 / 20 / 0.01453,
                80 => 0.287 / 20 / 0.01453,
                90 => 0.285 / 20 / 0.01453,
                )
    

end


@testset "halpin tsai" verbose = true begin
    Ef = 72.6
    Em = 1.5
    nu_f = 0.25
    nu_m = 0.39
    ar = 50
    vf = 0.1
    (Cht, constants) = SFRPMicroMechanics.halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar)
    display(Cht)
    display(constants)
    B = SFRPMicroMechanics.orientation_averaging_coefficients(Cht)
    Bref = (4680.46e-3, -18.18e-3, 16.97e-3, 2006.85e-3, 637.66e-3)

    for i in 1:5
        @info "B$i" B[i]
        @test B[i] ≈ Bref[i] atol=1e-4

    end

    el_const2 = SFRPMicroMechanics.extract_orthotropic_constants(Cht)
    display(el_const2)
    


end