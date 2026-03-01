@testset "Tensor Conversions" verbose = true begin
    C = SFRPMicroMechanics.isotropic_stiffness(210.0, 0.3)
    tens = SFRPMicroMechanics.convert_66_to_3333(C; mandel = false)
    Cc = SFRPMicroMechanics.convert_3333_to_66(tens;mandel = false)
    tensc = SFRPMicroMechanics.convert_66_to_3333(Cc; mandel = false)

    @test all(C .≈ Cc)
    @test all(tens .≈ tensc)


end

@testset "Basic Stiffness Conversion" verbose = true begin
    E, nu = 10.0, 0.3
    p =SFRPMicroMechanics.IsotropicElasticParameters(E, nu)
    C = SFRPMicroMechanics.stiffness_matrix_voigt(p)
    @test SFRPMicroMechanics.is_isotropic(C) #dooh
    ct = SFRPMicroMechanics.extract_orthotropic_constants(C)
    @test ct.E1 ≈ ct.E2 ≈ ct.E3 ≈ E
    @test ct.nu21  ≈ ct.nu32 ≈ ct.nu31 ≈ nu
    @test ct.G12 ≈ ct.G23 ≈ ct.G31 ≈ E/ (2 * (1 + nu))

    E1, E2, E3 = 20.0, 6.0, 4.0
    nu21 = 0.05 
    nu32 = 0.1
    nu31 = 0.25
    G12 = 8.0
    G23 = 6.4
    G31 = 4.5
    p = SFRPMicroMechanics.OrthotropicElasticParameters(;E1, E2, E3, nu21, nu32, nu31, G12, G23, G31) 
    C = SFRPMicroMechanics.stiffness_matrix_voigt(p)

    ct = SFRPMicroMechanics.extract_orthotropic_constants(C)

    @test ct.E1 ≈ E1
    @test ct.E2 ≈ E2
    @test ct.E3 ≈ E3
    @test ct.nu21 ≈ nu21
    @test ct.nu32 ≈ nu32
    @test ct.nu31 ≈ nu31
    @test ct.G12 ≈ G12
    @test ct.G23 ≈ G23
    @test ct.G31 ≈ G31

end



# @testset "SFRPMicroMechanics Physical Validation" verbose=true begin
@testset "No Fibers" verbose = true begin
    # Test 1: Zero fibers
    #vf = 0
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.0, 15.0, 0.7, 0.25
    C = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    props = SFRPMicroMechanics.extract_orthotropic_constants(C)
    # @info props
    @test isapprox(props.E1, Em, atol=1e-2)
    @test isapprox(props.nu21, num, atol=1e-4)
    @test isapprox(props.E2, Em, atol=1e-1)
    @test isapprox(props.G12, Em / (2*(1+num)), atol=1e-1)
end

@testset "Eshelby Tensor" verbose = true begin
    AR = 1e6 #S1111 -> 0 when AR-> Inf
    nu = 0.3
    shape =SFRPMicroMechanics.SpheroidalInclusion(nu, AR)
    S  = SFRPMicroMechanics.eshelby_tensor(shape) |> SFRPMicroMechanics.convert_3333_to_66
    # display(S)
    @test S[1,1] <= sqrt(eps(Float64))
    
    AR = 1e-10 #S1111 -> 1 when AR-> 0
    shape = SFRPMicroMechanics.SpheroidalInclusion(nu, AR) 
    S  = SFRPMicroMechanics.eshelby_tensor(shape) |> SFRPMicroMechanics.convert_3333_to_66
    @test S[1,1] ≈ 1
    Slimit = SFRPMicroMechanics.ThinDiscInclusion(nu, AR) |> SFRPMicroMechanics.eshelby_tensor |> SFRPMicroMechanics.convert_3333_to_66
    # @test all(Slimit .≈ S)
    @test Slimit[1,1] ≈ S[1,1]
end


@testset "Anisotropic Properties" verbose = true begin
    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.3, 15.0, 0.7, 0.25
    # C = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    # props = SFRPMicroMechanics.extract_orthotropic_constants(C)
    # @info props
    mandel = true
    # shape = SFRPMicroMechanics.SpheroidalInclusion(num, ar)
    Cavg = SFRPMicroMechanics.compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    props = SFRPMicroMechanics.extract_orthotropic_constants(Cavg)
    @test props.E1 > props.E2
    @test isapprox(props.E2, props.E3, rtol=1e-3)
    display(Cavg)
    
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
    mandel = true
    #mori tanaka should already be isotropic 
    pm = SFRPMicroMechanics.IsotropicElasticParameters(Em, num)
    Cm = SFRPMicroMechanics.stiffness_matrix_voigt(pm;mandel)
    pf = SFRPMicroMechanics.IsotropicElasticParameters(Ef, nuf)
    Cf = SFRPMicroMechanics.stiffness_matrix_voigt(pf;mandel)
    Cmt = SFRPMicroMechanics.mori_tanaka(Cm, Cf, vf, ar, num)
    @test SFRPMicroMechanics.is_structurally_isotropic(Cmt) #this is very stringent
    
    # display(Cmt)
    elprops = SFRPMicroMechanics.extract_orthotropic_constants(Cmt)
    @test elprops.E1 ≈ elprops.E2 ≈ elprops.E3
    @test elprops.nu21 ≈ elprops.nu32 ≈ elprops.nu31
    @test elprops.G12 ≈ elprops.G23 ≈ elprops.G31

    Ceff = compute_orthotropic_properties(Em, num, Ef, nuf, vf, ar, a11, a22)
    display(Ceff) #almost isotropic...
    elprops = SFRPMicroMechanics.extract_orthotropic_constants(Ceff)
    display(elprops)
    @test SFRPMicroMechanics.is_isotropic(Ceff)
end
    
@testset "Unidirectional (UD) Alignment" verbose = true begin
    # a11 = 1.0 means all fibers are perfectly parallel.
    # This must yield a Transversely Isotropic material.
    a11, a22  = 1.0,  0.0
    Em, num, Ef, nuf = 2500.0, 0.3, 70e3, 0.2
    aspect_ratio = 10.0
    vf = 0.17
    Ceff = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    res = props = SFRPMicroMechanics.extract_orthotropic_constants(Ceff)
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
    C = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    @test SFRPMicroMechanics.is_isotropic(C)
    
end

@testset "Bounds Check (Rule of Mixtures)" verbose = true begin
    # The longitudinal modulus E1 of a UD composite (a11=1) 
    # should never exceed the Voigt Upper Bound (Rule of Mixtures).
    a11, a22  = 1.0, 0.0
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 1000.0
    vf = 0.2
    C = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22) # AR=1000 ~ continuous
    res = SFRPMicroMechanics.extract_orthotropic_constants(C)
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
    @test SFRPMicroMechanics.apparent_modulus(p, angle) ≈ p.E1
    angle = 90.0
    @test SFRPMicroMechanics.apparent_modulus(p, angle) ≈ p.E2

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
    
    #3D results should be the same as 2D in the 1-2 plane
    angles = 0.0:15:360
    phi = 90
    for theta in angles
        E2d = SFRPMicroMechanics.apparent_modulus(p, theta)
        E3d = SFRPMicroMechanics.apparent_modulus(p, theta, phi)
        @test E2d ≈ E3d atol=1e-4
        # @info "E2d = $E2d,  E3d = $E3d"
    end

    C = SFRPMicroMechanics.stiffness_matrix_voigt(p)
    theta = rand()*360
    @test SFRPMicroMechanics.apparent_modulus(p, theta) ≈ SFRPMicroMechanics.apparent_modulus(C, theta)

end


@testset "halpin tsai" verbose = true begin
    """
        values from "A comparatice study between fiber orientation closure approximations and a new orthotropic closure"
        by Ahmad Al-Qudsi, Hakan Celik, Jonas Neuhaus, Christian Hopmann

    """
    Ef = 72.6
    Em = 1.5
    nu_f = 0.25
    nu_m = 0.39
    ar = 50
    vf = 0.1
    (Cht, constants) = SFRPMicroMechanics.halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar)
    
    display(constants)
    B = SFRPMicroMechanics.orientation_averaging_coefficients(Cht)
    Bref = (4680.46e-3, -18.18e-3, 16.97e-3, 2006.85e-3, 637.66e-3)

    for i in 1:5
        @info "B$i" B[i]
        @test B[i] ≈ Bref[i] atol=1e-4

    end

    el_const2 = SFRPMicroMechanics.extract_orthotropic_constants(Cht)
    display(el_const2)
    
    (a1, a2, a3, a4, a5, a6) = (0.7171, 0.2389, 0.0439, -0.0043, 0.0068, -0.0165)

    N2 =SFRPMicroMechanics.SMatrix{3,3}([a1 a6 a5;
                                        a6 a2 a4;
                                        a5 a4 a3])

    N4 = SFRPMicroMechanics.convert_3333_to_66(
                            SFRPMicroMechanics.HL2_closure(N2);mandel=true)
    # display(N4)

    # A11 = 0.6
    # A12 = 0.0955
    # A13 = 0.0217
    # A14 = -0.0015
    # A15 = 0.0049
    # A16 = -0.0084
    # A22 = 0.1343
    # A23 = 0.0091
    # A24 = -0.0019
    # A25 = 0.0009
    # A26 = -0.0072
    # A33 = 0.0131
    # A34 = -0.0009
    # A35 = 0.0009
    # A36 = -0.0009
    # A44 = 0.0091
    # A55  =0.0217
    # A66 = 0.0955
    # A45 = -0.0009
    # A46 = 0.0009
    # A56 = -0.0015

    # N4_ref = [A11 A12 A13 A14 A15 A16;
    #           A12 A22 A23 A24 A25 A26;
    #           A13 A23 A33 A34 A35 A36;
    #           A14 A24 A34 A44 A45 A46;
    #           A15 A25 A35 A45 A55 A56;
    #           A16 A26 A36 A46 A56 A66]

    # closure = SFRPMicroMechanics.hybrid_closure
    # C_avg = SFRPMicroMechanics.orientation_average(Cht, N2; closure)
    # @info "Halpin Tsai"
    # display(Cht)
    # @info "Averaged global Cs"
    # display(C_avg)
    # el_const3 = SFRPMicroMechanics.extract_orthotropic_constants(C_avg)
    # display(el_const3)

    # N2_mat = SFRPMicroMechanics.FullOrientationTensor(a1, a2, a4, a5, a6) |> 
    #          SFRPMicroMechanics.OrientationTensor |>
    #          SFRPMicroMechanics.to_matrix
    # C_avg_mat = SFRPMicroMechanics.orientation_average(Cht, N2_mat; closure)
    # display(C_avg_mat)
    # el_const4 = SFRPMicroMechanics.extract_orthotropic_constants(C_avg_mat)
    # display(el_const4)
    # for (a4ii, Aii) in zip(a4, A4_ref)
    #     @test a4ii ≈ Aii atol=1e-3
    # end

end


@testset "Thermal Expansion" verbose = true begin
    Em, num = 3.5, 0.35 
    Cm = SFRPMicroMechanics.IsotropicElasticParameters(Em ,num) |> 
         SFRPMicroMechanics.stiffness_matrix_voigt

    alfa_m = 100e-6
    cte_m = SFRPMicroMechanics.ThermalExpansion(alfa_m)

    Ef, nuf = 75.0, 0.22
    Cf = SFRPMicroMechanics.IsotropicElasticParameters(Ef, nuf) |>
         SFRPMicroMechanics.stiffness_matrix_voigt

    alfa_f = 5e-6
    cte_f = SFRPMicroMechanics.ThermalExpansion(alfa_f)
    AR = 50
    Sf = SFRPMicroMechanics.SpheroidalInclusion(num, AR) |>
         SFRPMicroMechanics.eshelby_tensor |>
         SFRPMicroMechanics.convert_3333_to_66

    vf = 0.2
    a11, a22 = 0.7, 0.2
    """
    matrix = (;stiffness = Cm, 
                thermal_expansion = αm)
    fiber_properties = (;stiffness = [Cfi],
                        thermal_expansion = [αfi],
                        eshelby_tensors = [Sfi],
                        volume_fractions = [vfi])
    """
    matrix = (;stiffness = Cm,
                thermal_expansion = cte_m
                )
    fibers = (;stiffness =[Cf],
                thermal_expansion = [cte_f],
                eshelby_tensors = [Sf],
                volume_fractions = [vf]
                )

    cte_eff = SFRPMicroMechanics.ThermalExpansion(matrix, fibers)
    display(cte_eff)
    @test alfa_f < cte_eff.alpha1 < alfa_m

    #might as well
    # ThermalExpansion(Em::T, num::T, alpham::T, Ef::T, nuf::T, alphaf::T, vf::T, AR::T, a11::T, a22::T)
    cte_eff = ThermalExpansion(Em, num, alfa_m, Ef, nuf, alfa_f, vf, AR, a11, a22)
    #dir 1 alfa is closer to fiber than dir 2 and 3
    display(cte_eff)
    @test abs(alfa_f - cte_eff.alpha1) < abs(alfa_f - cte_eff.alpha2)
end


@testset "Comparison to homopy" verbose = true begin
    #isotropic fibers
    mandel = true
    pf = SFRPMicroMechanics.IsotropicElasticParameters(242, 0.1)
    Cf = SFRPMicroMechanics.stiffness_matrix_voigt(pf; mandel)
    nu_m = 0.35
    pm = SFRPMicroMechanics.IsotropicElasticParameters(2.0, nu_m)
    Cm = SFRPMicroMechanics.stiffness_matrix_voigt(pm; mandel)
    ar = 20
    vf = 0.2
    Ceff = SFRPMicroMechanics.mori_tanaka(Cm, Cf, vf, ar, nu_m; 
                    fiber_shape = SpheroidalInclusion, 
                    mandel,
                    )
    # @info "Ceff"
    # display(Ceff)
    Ceff_homopy = [22.56165115 2.12961706 2.12961706 0 0 0;
                    2.12961706 4.28614074 2.20904525 0 0 0;
                    2.12961706 2.20904525 4.28614074 0 0 0;
                        0         0        0     2.0770955 0 0;
                        0         0        0         0  2.21728527 0;
                        0         0        0         0     0   2.21728527]
    # @info "Ceff homopy"
    # display(Ceff_homopy)
    @test all(isapprox.(Ceff, Ceff_homopy, atol = 1e-4))
    # for (c, c_hom) in zip(Ceff, Ceff_homopy)
    #     @test c ≈ c_hom atol=1e-4
    # end

    a11, a22 = 0.7, 0.25
    orientation_tensor =SFRPMicroMechanics.OrientationTensor(a11, a22) |> SFRPMicroMechanics.to_matrix
    Cavg = SFRPMicroMechanics.orientation_average(Ceff, orientation_tensor; mandel)

    @info "Cavg"
    display(Cavg)

    Cavg_homopy = [13.72300436  5.12826607    3.00247743    0       0         0;
                    5.12826607  5.80792549    2.42678122    0       0         0;
                    3.00257743  2.42678122    4.12680095    0       0         0;
                    0           0             0        2.25192461   0         0;
                    0           0             0             0   1.90998018    0;
                    0           0             0             0       0     3.20117434]
    @info "Cavg homopy"
    display(Cavg_homopy)
    @test all(isapprox.(Cavg, Cavg_homopy, atol =1e-4))
    # for (ca, ca_hom) in zip(Cavg, Cavg_homopy)
    #     @test ca ≈ ca_hom
    # end

    #test the hybrid closure element by element // also wrong
    # @info "N4 hybrid"
    # N4 = display(SFRPMicroMechanics.hybrid_closure(orientation_tensor))


    inclusion = SFRPMicroMechanics.SpheroidalInclusion(nu_m, ar)
    S = SFRPMicroMechanics.convert_3333_to_66( 
                        SFRPMicroMechanics.eshelby_tensor(inclusion) ;
                        mandel = true)
    # @info "S"
    # display(S)
    S_homopy = [1.52433537e-2 -6.13043062e-4 -6.13043062e-4       0           0           0;
                2.63166567e-1  6.90820632e-1  7.74657062e-2       0           0           0;
                2.63166567e-1  7.74657062e-2  6.90820632e-1       0           0           0;
                    0              0              0          6.13354926e-1    0           0;
                    0              0              0               0        4.94880228e-1  0;
                    0              0              0               0           0       4.9480228e-1]
    # @info "S_homopy"
    # display(S_homopy)
    @test all(isapprox.(S, S_homopy, atol=1e-4))
    # for (s, s_hom) in zip(S, S_homopy)
    #     @test s ≈ s_hom atol=1e-4
    # end

    #transverse isotropic fibers

    E1_c = 230.0
    E2_c = E3_c = 50.0
    G12_c = 10.0
    nu21_c = nu31_c = 0.03
    nu23_c = 0.39

    G23_c = E2_c / (2 * (1 +nu23_c))
    G13_c = G12_c

    pf = SFRPMicroMechanics.OrthotropicElasticParameters(;E1 = E1_c,
                                                        E2 = E2_c, E3 = E3_c,
                                                        G12 = G12_c, G23 = G23_c, G31 = G13_c,
                                                        nu21 = nu21_c, nu23 = nu23_c, nu31 = nu31_c)
    display(pf)
    Cf = SFRPMicroMechanics.stiffness_matrix_voigt(pf; mandel)
    # @info "Cf"
    # display(Cf)
    Cf_homopy = [233.16492721    11.46712757    11.46712757    0        0        0;
                  11.46712757    59.53317516    23.56195214    0        0        0;
                  11.46712757    23.56195214    59.53317516    0        0        0;
                   0             0               0          35.97122302 0        0;
                   0             0               0             0     20.0         0;
                   0             0               0             0        0       20.0]
    # @info "Cf_homopy"
    # display(Cf_homopy)
    @test all(isapprox.(Cf, Cf_homopy, atol = 1e-4))
    # for (c, c_hom) in zip(Cf, Cf_homopy)
    #     @test c ≈ c_hom
    # end

    Ceff = SFRPMicroMechanics.mori_tanaka(Cm, Cf, vf, ar, nu_m; 
                    fiber_shape = SpheroidalInclusion, 
                    mandel,
                    )
    @info "Ceff trans"
    display(Ceff)

    Ceff_homopy = [29.92384992  2.21785213    2.21785213    0     0    0 ;
                    2.21785213  4.54152923    2.32401817    0     0    0 ;
                    2.21785213  2.32401817    4.54152923    0     0    0 ;
                    0           0             0         2.21751107 0   0;
                    0           0             0             0  2.30146776  0;
                    0           0             0             0      0   2.30146776]

    @info "Ceff_homopy trans"
    display(Ceff_homopy)
    @test all(Ceff .≈ Ceff_homopy)
    # for (c, c_hom) in zip(Ceff, Ceff_homopy)
    #     @test c ≈ c_hom
    # end


end


@testset "Closures" verbose = true begin

    S = SFRPMicroMechanics
    CTs = [S.LinearClosure, 
           S.QuadraticClosure, 
           S.HybridClosure, 
           S.HL1Closure,
           S.HL2Closure,
           S.ORF, 
           S.ORL, 
           S.ORFM, 
           S.ORW, 
           S.ORW3,
           S.IBOF]

    Q, rot = qr(randn(3,3))
    R = Matrix(Q)
    adiag = diagm([0.6, 0.3, 0.1])
    amat = R * adiag * R'
    a = S.FullOrientationTensor(amat[1, 1], amat[2,2], amat[2,3], amat[1,3], amat[1,2])
    a2 = S.to_matrix(a)
    for CT in CTs
        @info "CT = $CT"
        # @test S.test_closure_approximation(a, CT; tol = 1e-8)
        N4 = S.closure(a, CT)
        @test S.is_normalized(N4,a2)
        @test S.is_symmetric(N4)
        # @test S.is_positive_semidef(N4)
    end


end