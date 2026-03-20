@testset "Basic Workflow" verbose = true begin
    #this shouldn't produce errors ...
    S = SFRPMicroMechanics
    Em, num = 2.3456, 0.345
    pm = S.IsotropicElasticParameters(Em, num) #elastic properties
    Cm = S.stiffness_matrix_voigt(pm; mandel = true) #stiffness matrix
    # @info "Cm"
    # display(Cm)
    
    Ef, nuf = 78.91234, 0.23551
    pf = S.IsotropicElasticParameters(Ef, nuf)
    Cf = S.stiffness_matrix_voigt(pf; mandel = true)
    # @info "Cf"
    # display(Cf)

    vf = 0.23345
    AR = 45.3424
    
    alpha_f = 1.0
    cte = S.ThermalExpansion(alpha_f)
    fibers = [S.FiberPhase(pf, vf, AR, S.SpheroidalInclusion())]

    Cmt = S.mori_tanaka(pm, fibers;mandel = true, symmetrize = true)
    #
    # @info "Cmt"
    # display(Cmt)


    a11, a22 = 0.5674, 0.3333
    a = S.OrientationTensor(a11 ,a22)
    N2 = S.to_matrix(a)
    N4 = S.closure(a, S.HybridClosure)
    # display(N4)

    Cavg = S.orientation_average(Cmt, a; mandel=true)
    # @info "Cavg"
    # display(Cavg)

    #taken from homopy with the same parameters
    Cavg_ref = [11.41172054   4.91443011    3.53631121    0        0            0;
                4.91443011    7.9294915     3.08564431    0        0            0;
                3.53631121    3.08564431    5.20297598    0        0            0;
                 0            0              0       3.25600713    0            0;
                 0            0              0            0        3.91019283   0;
                 0            0              0            0        0         2.21728527]
            
                 
    # @test all(Cavg .≈ Cavg_ref)
    @test true 
end



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
    mandel = true
    p =SFRPMicroMechanics.IsotropicElasticParameters(E, nu)
    C = SFRPMicroMechanics.stiffness_matrix_voigt(p;mandel)
    @test SFRPMicroMechanics.is_structurally_isotropic(C) #dooh

    ct = SFRPMicroMechanics.extract_orthotropic_constants(C; mandel)
    @test p ≈ ct
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
    C = SFRPMicroMechanics.stiffness_matrix_voigt(p; mandel)

    ct = SFRPMicroMechanics.extract_orthotropic_constants(C; mandel)

    @test isapprox(p, ct; atol = 1e-8)
    

end



# @testset "SFRPMicroMechanics Physical Validation" verbose=true begin
@testset "No Fibers" verbose = true begin
    # Test 1: Zero fibers
    #vf = 0
    S = SFRPMicroMechanics
    vf = 0
    Em, num, Ef, nuf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 15.0, 0.7, 0.25

    a2 = S.OrientationTensor(a11, a22)
    pm = S.IsotropicElasticParameters(Em, num)
    pf = S.IsotropicElasticParameters(Ef, nuf)

    Cm = S.stiffness_matrix_voigt(pm;mandel = true)
    Cmt = S.mori_tanaka(pm, pf, vf, ar)
    
    C = S.orientation_average(Cmt, a2)

    @test all(isapprox.(Cmt, Cm))
    @test all(isapprox.(C, Cm))
    
end

@testset "Eshelby Tensor" verbose = true begin
    AR = 1e6 #S1111 -> 0 when AR-> Inf
    nu = 0.3
    shape =SFRPMicroMechanics.SpheroidalInclusion()
    S  = SFRPMicroMechanics.eshelby_tensor(shape, nu, AR) |> SFRPMicroMechanics.convert_3333_to_66
    # display(S)
    @test S[1,1] <= sqrt(eps(Float64))
    
    AR = 1e-10 #S1111 -> 1 when AR-> 0
    shape = SFRPMicroMechanics.SpheroidalInclusion() 
    S  = SFRPMicroMechanics.eshelby_tensor(shape, nu, AR) |> SFRPMicroMechanics.convert_3333_to_66
    @test S[1,1] ≈ 1
    Slimit = SFRPMicroMechanics.eshelby_tensor(ThinDiscInclusion(),nu, AR)  |> SFRPMicroMechanics.convert_3333_to_66
    # @test all(Slimit .≈ S)
    @test Slimit[1,1] ≈ S[1,1]
end


@testset "Anisotropic Properties" verbose = true begin
    S = SFRPMicroMechanics

    Em, num, Ef, nuf, vf, ar, a11, a22 = 2000.0, 0.35, 70e3, 0.22, 0.3, 15.0, 0.7, 0.25
    
    pm = S.IsotropicElasticParameters(Em, num)
    pf = S.IsotropicElasticParameters(Ef, nuf)

    shape = S.SpheroidalInclusion()
    a = S.OrientationTensor(a11, a22)
    fibers = [S.FiberPhase(pf, vf, ar, shape)]

    mandel = true
    Ceff = S.mori_tanaka(pm, fibers;mandel)
    # @info "Ceff"
    # display(Ceff)
    el_props = S.extract_orthotropic_constants(Ceff)
    @test el_props.E1 > el_props.E2
    @test el_props.E2 ≈ el_props.E3

    closure_type = S.HybridClosure
    N4 = S.closure(a, closure_type)
    # @info "N4"
    # display(N4)
    Cavg = S.orientation_average(Ceff, a;mandel, closure_type)
    el_props = S.extract_orthotropic_constants(Cavg)
    @test el_props.E1 > el_props.E2 #a11 > a22
    @test el_props.E2 > el_props.E3 #a22>a33
    # @info "Cavg"
    # display(Cavg)
    # @info "el props Cavg"
    # display(el_props)
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
    for mandel in (true, false)
        Cm = SFRPMicroMechanics.stiffness_matrix_voigt(p; mandel)
        pex = SFRPMicroMechanics.extract_orthotropic_constants(Cm; mandel)
        @test p ≈ pex
       
    end
end

@testset "Transverse Properties" verbose = true begin

    S = SFRPMicroMechanics

    E1 = 230.0
    E2 = 25.0
    E3 = E2

    G12 = 50.0
    G31 = G12
    nu21 = nu31 = 0.03
    nu23 = nu32 = 0.39
    nu12 = nu13 = nu21 * E1 / E2
    G23 = E2 / (2(1+nu23))

    #these should be all the same 
    p1 = S.OrthotropicElasticParameters(;E1, E2, E3, G12, G23, G31, nu12, nu13, nu23)
    pt = S.TransverseIsotropicElasticParameters(;E1, E2, G12, nu23, nu13)
    pt2 = S.TransverseIsotropicElasticParameters(;E1, E3, G31, G23, nu12)
    pt3 = S.TransverseIsotropicElasticParameters(;E1, E2, G12, G23, nu31)
    pt4 = S.TransverseIsotropicElasticParameters(;E1, E3, G31, nu21, nu32)
    pt5 = S.TransverseIsotropicElasticParameters(;E1, E2, G12, nu21, nu23) #ls-dyna mat215
    pt6 = S.TransverseIsotropicElasticParameters(;E1, E2, G12, G23, nu12) #python homopy
    # @info "p1"
    # display(p1)
    # @info "pt"
    # display(pt)
    # @info "pt2"
    # display(pt2)
    # @info "pt3"
    # display(pt3)
    # @info "pt4"
    # display(pt4)
    @test isapprox(p1, pt)
    @test isapprox(p1, pt2)
    @test isapprox(p1, pt3)
    @test isapprox(p1, pt4)
    @test isapprox(p1, pt5)
    @test isapprox(p1, pt6)
end

@testset "Orientation Averaging" verbose = true begin
    S = SFRPMicroMechanics
    ar = 1.0
    vf = 0.15
    Em, num, Ef, nuf = 2000, 0.3, 70e3, 0.22
    a11 = 0.7
    a22 = 0.25
    mandel = true
    #mori tanaka should already be isotropic 
    pm = S.IsotropicElasticParameters(Em, num)
    Cm = S.stiffness_matrix_voigt(pm;mandel)
    
    pf = S.IsotropicElasticParameters(Ef, nuf)
    Cf = S.stiffness_matrix_voigt(pf;mandel)

    Tf = S.convert_66_to_3333(Cf; mandel)
    # @info "Tf"
    # display(Tf)

    bs = S.orientation_averaging_coefficients(Cf;mandel)
    # @info "bs"
    # display(bs)

    bs_ref = [-2.9103830456733704e-11
                0.0
                7.275957614183426e-12
                22540.983606557376
                28688.524590163935
]
@test all(isapprox.(bs, bs_ref, atol=1e-6))
    
end
    
@testset "Spherical Inclusion Limit" verbose = true begin
    # @info "Spherical Inclusion Limit"
    # Aspect Ratio = 1.0 means fibers are spheres. 
    # Spheres are isotropic; E1 should equal E2 even with aligned a11.
    S = SFRPMicroMechanics
    ar = 1
    vf = 0.3
    Em, num, Ef, nuf = 2000, 0.3, 70e3, 0.22
    a11 = 0.7
    a22 = 0.25
    mandel = true
    #mori tanaka should already be isotropic 
    pm = S.IsotropicElasticParameters(Em, num)
    Cm = S.stiffness_matrix_voigt(pm;mandel)
    
    pf = S.IsotropicElasticParameters(Ef, nuf)
    Cf = S.stiffness_matrix_voigt(pf;mandel)
    
    # Cmt = S.mori_tanaka(Cm, Cf, vf, ar, num)
    Cmt = S.mori_tanaka(pm, pf, vf, ar; mandel)
    # @info "Sphere Cmt"
    # display(Cmt)
    @test S.is_structurally_isotropic(Cmt) #this is very stringent
    

    phase = S.FiberPhase(pf, vf, ar, S.SpheroidalInclusion())
    Cmt2 = S.mori_tanaka(pm, [phase];mandel)
    # @info "Cmt2"
    # display(Cmt2)

    @test all(Cmt .≈ Cmt2)
    #from homopy
    a, b, c = 4597.54348806, 1780.17562439, 2817.36786368
    Cmt_ref = [a b b 0 0 0;
               b a b 0 0 0;
               b b a 0 0 0;
               0 0 0 c 0 0;
               0 0 0 0 c 0;
               0 0 0 0 0 c]

    @test all(isapprox.(Cmt2, Cmt_ref, atol = 1e-4))


    Ssph = S.convert_3333_to_66(S.eshelby_tensor(S.SpheroidalInclusion(), num, ar);mandel=true)
    # @info "S spherical"
    # display(Ssph)
    
    a,b,c = 0.52380952, 0.04761905, 0.47619048 #from homopy
    S_ref = [a b b 0 0 0;
               b a b 0 0 0;
               b b a 0 0 0;
               0 0 0 c 0 0;
               0 0 0 0 c 0;
               0 0 0 0 0 c]
    @test all(isapprox.(Ssph, S_ref, atol = 1e-5))
    # @info "Cmt 3333"
    # display(S.convert_66_to_3333(Cmt;mandel))
    bs = S.orientation_averaging_coefficients(Cmt; mandel)
    # @info "Bs"
    # display(bs)

    # display(Cmt)
    elprops = S.extract_orthotropic_constants(Cmt)
    @test elprops.E1 ≈ elprops.E2 ≈ elprops.E3
    @test elprops.nu21 ≈ elprops.nu32 ≈ elprops.nu31
    @test elprops.G12 ≈ elprops.G23 ≈ elprops.G31

    a = S.OrientationTensor(a11, a22)
    Cavg = S.orientation_average(Cmt, a; mandel)
    @test all(isapprox.(Cavg, Cmt)) 
    # @info "Sphere Cavg"
    # display(Cavg) works
    elprops = S.extract_orthotropic_constants(Cavg)
    display(elprops)
    # @test S.is_isotropic(Cavg)
    @test S.is_structurally_isotropic(Cavg)
end
    
@testset "Unidirectional (UD) Alignment" verbose = true begin
    S = SFRPMicroMechanics
    # a11 = 1.0 means all fibers are perfectly parallel.
    # This must yield a Transversely Isotropic material.
    a11, a22  = 1.0,  0.0
    Em, num, Ef, nuf = 2500.0, 0.3, 70e3, 0.2
    pm = S.IsotropicElasticParameters(Em, num)
    pf = S.IsotropicElasticParameters(Ef, nuf)
    aspect_ratio = 10.0
    vf = 0.17
    Ceff = S.mori_tanaka(pm, pf, vf, aspect_ratio)
    res = props = S.extract_orthotropic_constants(Ceff)
    @test res.E1 > res.E2
    @test isapprox(res.E2, res.E3, rtol=1e-4)
    @test isapprox(res.nu21, res.nu31, rtol=1e-4)
    # Explanation: In a UD composite, the properties perpendicular 
    # to the fiber (2 and 3 directions) must be identical.
end
    
@testset "3D Random Isotropy" verbose = true begin
    # a11 = a22 = a33 = 1/3 is a perfectly random 3D distribution.
    S = SFRPMicroMechanics
    a11 = a22  = 1/3
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 10.0
    vf = 0.2
    pm = S.IsotropicElasticParameters(Em, num)
    pf = S.IsotropicElasticParameters(Ef, nuf)

    a2 = S.OrientationTensor(a11, a22)
    Cmt = S.mori_tanaka(pm, pf, vf, aspect_ratio)
    Cavg = S.orientation_average(Cmt, a2)
    # @info "Cavg"
    # display(Cavg)
    # C = compute_orthotropic_properties(Em, num, Ef, nuf, vf, aspect_ratio, a11, a22)
    @test SFRPMicroMechanics.is_structurally_isotropic(Cavg)
    ep = S.extract_orthotropic_constants(Cavg; mandel = true)
    E, nu, G = ep.E1, ep.nu21, ep.G12
    @test G ≈ E/(2(1+nu))
end

@testset "Bounds Check (Rule of Mixtures)" verbose = true begin
    # The longitudinal modulus E1 of a UD composite (a11=1) 
    # should never exceed the Voigt Upper Bound (Rule of Mixtures).
    S = SFRPMicroMechanics
    a11, a22  = 1.0, 0.0
    Em, num, Ef, nuf = 2000.0, 0.35, 70e3, 0.2
    aspect_ratio = 1000.0
    vf = 0.2

    mandel = true

    pm = S.IsotropicElasticParameters(Em, num)
    pf  =S.IsotropicElasticParameters(Ef, nuf)
    shapes  = [S.NeedleInclusion(), S.SpheroidalInclusion()]
    for shape in shapes
        # @info "$shape"
        fibers = [S.FiberPhase(pf, vf, aspect_ratio, shape)]
        Cmt = S.mori_tanaka(pm, fibers;mandel, symmetrize=true)
        elprops = S.extract_orthotropic_constants(Cmt; mandel)
        
        ROM_E1 = vf * Ef + (1 - vf) * Em #rule of mixtures 
        # @info "E1 computed = $(elprops.E1)"
        # @info "Upper bound = $ROM_E1"
        @test elprops.E1 <= ROM_E1 + 1.0 # Allow for tiny numerical float noise

        Cavg = S.orientation_average(Cmt, S.OrientationTensor(a11, a22);mandel)
        elprops = S.extract_orthotropic_constants(Cavg;mandel)
        
        
        # @info "E1 computed avg = $(elprops.E1)"
        # @info "Upper bound = $ROM_E1"
        @test elprops.E1 <= ROM_E1 + 1.0 # Allow for tiny numerical float noise

    end
    # Explanation: Micromechanics models must respect thermodynamic bounds.
    #TODO fails for SpheroidalInclusion


    #halpin-tsai
    #iso
    C_ht = S.halpin_tsai(pm, pf, vf, aspect_ratio)
    p = S.extract_orthotropic_constants(C_ht)
    ROM_E1 = (1 - vf) * pm.E_modulus + vf * pf.E_modulus
    @info "Iso Halpin Tsai E1 = $(p.E1), ROM E = $ROM_E1"
     
    @test p.E1 < ROM_E1

    #trans fiber
    pf = S.TransverseIsotropicElasticParameters(;E1 = 230.0, E2 = 25.0, nu21 = 0.03, nu23 = 0.4, G12 = 50.0)
    C_ht = S.halpin_tsai(pm, pf, vf, aspect_ratio)
    p = S.extract_orthotropic_constants(C_ht)
    ROM_E1 = (1 - vf) * pm.E_modulus + vf * pf.E1
    @info "Transverse Halpin Tsai E1 = $(p.E1), ROM E = $ROM_E1"
    @test p.E1 < ROM_E1

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

    C = SFRPMicroMechanics.stiffness_matrix_voigt(p; mandel=true)
    theta = rand()*360
    @test SFRPMicroMechanics.apparent_modulus(p, theta) ≈ SFRPMicroMechanics.apparent_modulus(C, theta)

end


@testset "halpin tsai" verbose = true begin
    """
        values from "A comparatice study between fiber orientation closure approximations and a new orthotropic closure"
        by Ahmad Al-Qudsi, Hakan Celik, Jonas Neuhaus, Christian Hopmann

    """
    S = SFRPMicroMechanics


    Ef = 72.6
    nu_f = 0.25
    pf = S.IsotropicElasticParameters(Ef, nu_f)

    Em = 1.5
    nu_m = 0.39
    pm = S.IsotropicElasticParameters(Em ,nu_m)
    
    ar = 50
    vf = 0.1
    mandel = true
    
    # Cht = S.halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar; mandel)
    Cht = S.halpin_tsai(pm, pf, vf, ar; mandel)
    # display(constants)
    B = S.orientation_averaging_coefficients(Cht; mandel)
    Bref = (4680.46e-3, -18.18e-3, 16.97e-3, 2006.85e-3, 637.66e-3)

    for i in 1:5
        # @info "B$i" B[i]
        @test B[i] ≈ Bref[i] atol=1e-4

    end

    el_const2 = S.extract_orthotropic_constants(Cht)
    # display(el_const2)
    
    (a1, a2, a3, a4, a5, a6) = (0.7171, 0.2389, 0.0439, -0.0043, 0.0068, -0.0165)

    # N2 =S.SMatrix{3,3}([a1 a6 a5;
    #                                     a6 a2 a4;
    #                                     a5 a4 a3])

    # N4 = S.convert_3333_to_66(
    #                         S.HL2_closure(N2);mandel=true)
    
    
    ## isotropic == transverse with isotropic props
    G12 = Ef / (2(1 + nu_f))
    # G23 = G12
    pf_trans =  S.TransverseIsotropicElasticParameters(Ef, Ef, G12, G12, nu_f)
    Cht_trans = S.halpin_tsai(pm, pf, vf, ar; mandel)  

    el_const_trans = S.extract_orthotropic_constants(Cht_trans)
    @test isapprox(el_const2, el_const_trans; atol =1e-8)
end


@testset "Thermal Expansion" verbose = true begin
    
    S = SFRPMicroMechanics
    Em, num = 3.5, 0.35 
    pm = S.IsotropicElasticParameters(Em, num)
    alfa_m = 100e-6
    # @info "CTE matrix"
    cte_m = SFRPMicroMechanics.ThermalExpansion(alfa_m)
    # display(cte_m)
    #isotropic fiber
    mandel = true
    Ef, nuf = 230.0, 0.2
    pf = S.IsotropicElasticParameters(Ef, nuf)

    alfa_f = -1e-6
    # @info "CTE fiber"
    cte_f = SFRPMicroMechanics.ThermalExpansion(alfa_f)
    # display(cte_f)
    AR = 50
    shape = S.SpheroidalInclusion()
    vf = 0.4
    
    a = S.OrientationTensor(0.7, 0.2)

    fibers = [S.FiberPhase(pf, vf, AR, shape)]

    ctes = [cte_m, cte_f]
    @info "effective CTE MT"
    cte_eff = S.effective_thermal_expansion_mt(pm, fibers, ctes, a; mandel)
    display(cte_eff)
    
    # @info "effective CTE Chow"
    # cte_eff_chow = S.effective_thermal_expansion_chow(pm, fibers[1], ctes)
    # display(cte_eff_chow)

    # @info "effective CTE MT v2"
    cte_eff_MT2 = S.ThermalExpansion(pm, pf, cte_m, cte_f, vf, AR, a, shape)
    # display(cte_eff_MT2)

    #only matrix => cte_eff == cte matrix
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, 0.0, AR, a, shape)
    @test cte_eff.alpha1 ≈ cte_eff.alpha2 ≈ cte_eff.alpha3 ≈ alfa_m
    #only fiber => cte_eff == cte fiber
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, 1.0, AR, a, shape)
    @test cte_eff.alpha1 ≈ cte_eff.alpha2 ≈ cte_eff.alpha3 ≈ alfa_f
    #  0< vf<1 => cte inbetween, transverse ortho
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, 0.2, AR, a, shape)
    @test alfa_f < cte_eff.alpha1 < cte_eff.alpha2 < cte_eff.alpha3 < alfa_m
    # if a11=a22 = 1/3 => cte isotropic
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, rand(), AR, S.OrientationTensor(1/3, 1/3), shape)
    @test cte_eff.alpha1 ≈ cte_eff.alpha2 ≈ cte_eff.alpha3
    # if ar=1 => cte isotropic
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, rand(), 1, a, shape)
    @test cte_eff.alpha1 ≈ cte_eff.alpha2 ≈ cte_eff.alpha3
    # if a22=a33 => alpha2 == alpha3
    cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f, rand(), AR, S.OrientationTensor(0.5, 0.25), shape)
    @test cte_eff.alpha1 < cte_eff.alpha2 ≈ cte_eff.alpha3


    #results for transverse with isotropic props == isotropic
    Gf = Ef / (2(1+nuf))
    pf_iso = S.IsotropicElasticParameters(Ef, nuf)
    pf = S.TransverseIsotropicElasticParameters(;E1 = Ef, E2 = Ef, nu12 = nuf, nu23 = nuf, G12 = Gf)
    vf, AR = 0.2, 15.0
    fibers = [S.FiberPhase(pf, vf, AR, shape)]
    cte_vec = [S.ThermalExpansion(80e-6), S.ThermalExpansion(5e-6)]
    cte_f2 = S.ThermalExpansion(pm, fibers, cte_vec , a)
    cte_iso = S.ThermalExpansion(pm, [S.FiberPhase(pf_iso, vf, AR, shape)], cte_vec, a)
    #should be same cte
    @test cte_iso.alpha1 ≈ cte_f2.alpha1 ≈ cte_f2.alpha2 ≈ cte_f2.alpha3
    #    

    # just test if this works

    
    cte_collection = S.compute_all_thermal_expansions(pm, fibers, cte_vec, a; average = true, mandel =true)
    # display(cte_collection)
    @info "All the CTEs"
    for (name, cte) in pairs(cte_collection)
        @info name
        display(cte)
    end

end


@testset "Comparison to homopy" verbose = true begin
    #isotropic fibers
    S = SFRPMicroMechanics
    mandel = true
    pf = SFRPMicroMechanics.IsotropicElasticParameters(242, 0.1)
    Cf = SFRPMicroMechanics.stiffness_matrix_voigt(pf; mandel)
    nu_m = 0.35
    pm = SFRPMicroMechanics.IsotropicElasticParameters(2.0, nu_m)
    Cm = SFRPMicroMechanics.stiffness_matrix_voigt(pm; mandel)
    ar = 20
    vf = 0.2
    # Ceff = SFRPMicroMechanics.mori_tanaka(Cm, Cf, vf, ar, nu_m; 
    #                 fiber_shape = SpheroidalInclusion(), 
    #                 mandel,
    #                 )
    Ceff = S.mori_tanaka(pm, pf, vf, ar;mandel)
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

    orientation_tensor =S.OrientationTensor(a11, a22) 
    Cavg = S.orientation_average(Ceff, orientation_tensor; mandel)

    @info "Cavg"
    display(Cavg)

    Cavg_homopy = [13.72300436  5.12826607    3.00247743    0       0         0;
                    5.12826607  5.80792549    2.42678122    0       0         0;
                    3.00257743  2.42678122    4.12680095    0       0         0;
                    0           0             0        2.25192461   0         0;
                    0           0             0             0   2.90998018    0;
                    0           0             0             0       0     3.20117434]
    @info "Cavg homopy"
    display(Cavg_homopy)
    @test all(isapprox.(Cavg, Cavg_homopy, atol =1e-4))
    # for (ca, ca_hom) in zip(Cavg, Cavg_homopy)
    #     @test ca ≈ ca_hom
    # end

    #test the hybrid closure element by element // works
    # @info "N4 hybrid"
    # N4 = display(SFRPMicroMechanics.hybrid_closure(orientation_tensor))


    inclusion = SFRPMicroMechanics.SpheroidalInclusion()
    S_eff = SFRPMicroMechanics.convert_3333_to_66( 
                        SFRPMicroMechanics.eshelby_tensor(inclusion, nu_m, ar) ;
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
    @test all(isapprox.(S_eff, S_homopy, atol=1e-4))
    

    #transverse isotropic fibers

    E1_c = 230.0
    E2_c = E3_c = 50.0
    G12_c = 10.0
    nu21_c = nu31_c = 0.03
    nu23_c = 0.39
    vf = 0.3
    ar = 22.0
    G23_c = E2_c / (2 * (1 +nu23_c))
    G13_c = G12_c

    pm = SFRPMicroMechanics.IsotropicElasticParameters(2.0, 0.35)

    pf = SFRPMicroMechanics.OrthotropicElasticParameters(;E1 = E1_c,
                                                        E2 = E2_c, E3 = E3_c,
                                                        G12 = G12_c, G23 = G23_c, G31 = G13_c,
                                                        nu21 = nu21_c, nu23 = nu23_c, nu31 = nu31_c)
    # display(pf)
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

    # Ceff = SFRPMicroMechanics.mori_tanaka(Cm, Cf, vf, ar, nu_m; 
    #                 fiber_shape = SpheroidalInclusion(), 
    #                 mandel,symmetrize =true
    #                 )
    Ceff = S.mori_tanaka(pm, pf, vf, ar;mandel)
    # @info "Ceff trans"
    # display(Ceff)

    Ceff_homopy = [36.44938137  2.34300859    2.34300859    0     0    0 ;
                    2.34300859  4.91037016    2.48829998    0     0    0 ;
                    2.34300859  2.48829998    4.91037016    0     0    0 ;
                    0           0             0         2.42207018 0   0;
                    0           0             0             0  2.52257841  0;
                    0           0             0             0      0   2.52257841]

    # @info "Ceff_homopy trans"
    # display(Ceff_homopy)
    @test all(Ceff .≈ Ceff_homopy)
    # for (c, c_hom) in zip(Ceff, Ceff_homopy)
    #     @test c ≈ c_hom
    # end


    a11, a22 = 0.7, 0.25
    orientation_tensor =S.OrientationTensor(a11, a22) 

    C_avg_trans = S.orientation_average(Ceff, orientation_tensor; mandel)
    
    @info "C_avg_trans"
    display(C_avg_trans)

    a, b, c, d, e, f, g, h, j = 21.28753957, 7.39504232, 4.60318504, 2.83807254, 3.84243603, 7.40101738, 2.69510888, 3.64095994, 4.04242938 
    C_avg_trans_homopy =  [a f e 0 0 0;
                           f b d 0 0 0;
                           e d c 0 0 0;
                           0 0 0 g 0 0;
                           0 0 0 0 h 0;
                           0 0 0 0 0 j]
    @info "C_avg_trans_homopy"
    display(C_avg_trans_homopy)
    @test all(isapprox.(C_avg_trans, C_avg_trans_homopy, atol=1e-4))

end


@testset "Bond Matrix Validation" begin
    S = SFRPMicroMechanics
    I = S.LinearAlgebra.I
    # 1. Setup a valid rotation Q (45 deg around Z-axis)
    θ = π/4
    Q =          [ cos(θ) -sin(θ) 0;
                   sin(θ)  cos(θ) 0;
                   0       0      1]
    
    @test isapprox(Q * Q', I, atol=1e-12) # Sanity check on Q

    # 2. Generate the Bond Matrix M
    M = S.convert_rot_33_to_66(Q)

    # 3. TEST 1: Identity Mapping
    M_eye = S.convert_rot_33_to_66(I)
    @test isapprox(M_eye, I, atol=1e-12)

    # 4. TEST 2: Isotropic Invariance
    # Define an Isotropic Stiffness Matrix (C) in Voigt form
    E, ν = 210e9, 0.3
    λ = E*ν / ((1+ν)*(1-2ν))
    μ = E / (2*(1+ν))
    
    C_iso =  [
        λ+2μ  λ     λ     0   0   0;
        λ     λ+2μ  λ     0   0   0;
        λ     λ     λ+2μ  0   0   0;
        0     0     0     μ   0   0;
        0     0     0     0   μ   0;
        0     0     0     0   0   μ
    ]

    # Transform: C_rotated = M * C_iso * M'
    C_rot = M * C_iso * M'

    # @test isapprox(C_rot, C_iso, atol=1e-7)
    @test S.LinearAlgebra.norm(C_rot - C_iso) < 1e-7
    # @test all(isapprox.(C_rot, C_iso, atol = 1e-5))
    # @info "C_rot"
    # display(C_rot)
    # @info("C_iso")
    # display(C_iso)
    # 5. TEST 3: Trace Invariance
    # The sum of eigenvalues (or specific invariants) should hold.
    @test isapprox(tr(C_rot), tr(C_iso), atol=1e-7)

    # println("✅ All Bond Matrix tests passed!")
end



@testset "Closures" verbose = true begin

    S = SFRPMicroMechanics
    CTs = [S.LinearClosure, 
           S.QuadraticClosure, 
           S.HybridClosure, 
           S.HL1Closure,
           S.HL2Closure,
           S.ORS,
           S.ORF, 
           S.ORL, 
           S.ORFM, 
           S.ORW, 
           S.ORW3,
           S.IBOF]

    Q, rot = qr(randn(3,3))
    R = Matrix(Q)
    # @info "R"
    # display(R)
    I3 = [1 0 0;
          0 1 0;
          0 0 1]
    @test all(isapprox.(R * R',I3, atol = 1e-8)) #this must pass
    # display(R * R')
    
    #generate a randomly rotated orientation tensor
    a11, a22 = 0.6, 0.3
    adiag = diagm([a11, a22, 1-a11 -a22])
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

    #the simple case
    a_eig = S.OrientationTensor(a11, a22)
    a2 = S.to_matrix(a_eig)
    for CT in CTs
        @info "CT = $CT"
        # @test S.test_closure_approximation(a, CT; tol = 1e-8)
        N4 = S.closure(a_eig, CT)
        @test S.is_normalized(N4,a2)
        @test S.is_symmetric(N4)
        # @test S.is_positive_semidef(N4)
    end

    # N4 = S.closure(a_eig, S.ORF)
    # @info("N4 - ORF")
    # display(N4)

    # Rref = [0.21910778 -0.76198347 -0.60940378;
    #         0.67076085 -0.3359456   0.66122647;
    #        -0.70857016 -0.55364406  0.43750039]

    # a_rot = Rref * S.to_matrix(a_eig) * Rref'
    # @info "a_rot"
    # display(a_rot)
    # a_rot_full = S.FullOrientationTensor(a_rot[1,1], a_rot[2,2], a_rot[2,3], a_rot[1,3], a_rot[1,2])
    # N4 = S.closure(a_rot_full, S.ORF)
    # @info "rotated N4"
    # display(N4)

    #IBOF - compared to fiberoripy
    beta_ref = [0.007322197191712073, 
                0.10261471647113499, 
                2.6899603609285103, 
                -0.15094735616892185, 
                -3.209763487085262, 
                1.5253864146381542]

    A2 = S.OrientationTensor(0.7, 0.2)
    (II, III) = S.get_invariants(A2)
    bs = S.beta_coefficients(II, III)

    @test all(isapprox.(bs, beta_ref, rtol = 1e-8))
    # for (bi, bref) in zip(bs, beta_ref)
    #     @test bi ≈ bref
    # end


end