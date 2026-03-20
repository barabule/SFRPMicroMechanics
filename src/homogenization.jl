
"""
    FiberPhase - collects fiber properties:
        - elastic_properties - (isotropic or Orthotropic)
        - volume fraction of this phase
        - aspect ratio (ratiio of longest to shortest dim)
        - shape : InclusionGeometry (one of :SphericalInclusion(), SpheroidalInclusion(), NeedleInclusion(),
                    DiscInclusion(), ThinDiscInclusion()

"""
Base.@kwdef struct FiberPhase{EP<:AbstractElasticParameters, T<:Real, IG<:InclusionGeometry}   
    elastic_properties::EP
    volume_fraction::T
    aspect_ratio::T
    shape::IG

    function FiberPhase(ep, vf, ar, shap)
        VF, AR = promote(vf, ar)
        T = typeof(VF)
        EP = typeof(ep)
        IG = typeof(shap)
        @assert IG<:InclusionGeometry
        @assert EP<:AbstractElasticParameters
        return new{EP, T, IG}(ep, VF, AR, shap)
    end
end




"""
    mori_tanaka(pm::IsotropicElasticParameters, fibers::AbstractVector{<:FiberPhase}; 
                    mandel = false,
                    symmetrize = false)

Inputs:
    pm - elastic parameters of the matrix - IsotropicElasticParameters(E, nu)
    fibers - list of fibers [FiberPhase]

    mandel - stiffness matrices are Mandel form
    symmetrize - if true, symmetrizes the final stiffness matrix

Returns an effective stiffness matrix 6x6
"""
function mori_tanaka(pm::IsotropicElasticParameters, fibers::AbstractVector{<:FiberPhase}; 
                    mandel = true,
                    symmetrize = false)

    Cm = stiffness_matrix_voigt(pm; mandel)
    invCm = inv(Cm)
    I6 = I(6) # Identity matrix for Voigt/Mandel space

    # Calculate Matrix volume fraction
    vm = 1.0 - sum(f.volume_fraction for f in fibers)
    
    # Initialize terms for the summation
    # We'll accumulate X and Y directly to avoid storing large vectors
    X = zeros(eltype(Cm), 6, 6)
    ΣvfTr = zeros(eltype(Cm), 6, 6)
    
    nu_m = pm.nu

    for fiber in fibers
        vf = fiber.volume_fraction
        shape = fiber.shape
        AR = fiber.aspect_ratio
        
        # 1. Get Eshelby Tensor and transform to Voigt/Mandel
        S_tensor = eshelby_tensor(shape, nu_m, AR)
        Sr = convert_3333_to_66(S_tensor; mandel)
        
        # 2. Get Fiber Stiffness
        Cf = stiffness_matrix_voigt(fiber.elastic_properties; mandel)
        
        # 3. Calculate the Concentration Tensor (Tr)
        # Tr = [I + Sr * inv(Cm) * (Cf - Cm)]⁻¹
        Tr = inv(I6 + Sr * invCm * (Cf - Cm))
        
        # 4. Accumulate terms
        ΣvfTr += vf * Tr
        X += vf * (Cf - Cm) * Tr
    end
    
    # Y = vm*I + Σ(vf_r * Tr)
    Y = (vm * I6 + ΣvfTr)
    
    # Effective Stiffness: C_MT = Cm + X * Y⁻¹
    C_MT = Cm + X * inv(Y)

    if symmetrize
        C_MT = 0.5 * (C_MT + C_MT')
    end

    return C_MT
end

function mori_tanaka(pm::IsotropicElasticParameters, pf::AbstractElasticParameters, vf::Real, ar::Real;
                    shape = SpheroidalInclusion(), 
                    mandel = true,
                    symmetrize = false)

    fiber_phase = FiberPhase(pf, vf, ar, shape)

    return mori_tanaka(pm, [fiber_phase]; mandel, symmetrize)
end




function halpin_tsai(pm::IsotropicElasticParameters, pf::IsotropicElasticParameters,
                    volume_fraction, aspect_ratio; mandel = true)

    Em, nu_m = pm.E_modulus, pm.nu
    Ef, nu_f = pf.E_modulus, pf.nu
    vf = volume_fraction
    ar = aspect_ratio
    return halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar; mandel)
end

function halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar; mandel = true)

    ζ1 = 2ar
    η1 = (Ef/Em - 1)/ (Ef/Em + ζ1)

    E11_UD = Em * (1 + ζ1 * η1 * vf) / (1 - η1 * vf)
    # @info "E11 $E11_UD"
    nu12_UD = nu_f * vf + nu_m * (1 - vf)

    ζ2 = 2
    η2 = (Ef/Em - 1) / (Ef/Em + ζ2)

    E22_UD = E33_UD = Em * (1 + ζ2 * η2 * vf) / (1 - η2 * vf)
    # @info "E22 $E22_UD"
    ζ3 = 1
    Gf = Ef / (2 * (1 + nu_f))
    Gm = Em / (2 * (1 + nu_m))
    η3 = (Gf / Gm - 1) / (Gf / Gm + ζ3)
    
    G31_UD = G12_UD = Gm * (1 + ζ3 * η3 * vf) / (1 - η3 * vf)


    ζ4 = (1 + nu_m) / (3 - nu_m - 4nu_m^2)
    η4 = (Gf/Gm - 1) / (Gf/Gm + ζ4)
    
    G23_UD = Gm * (1 + ζ4 * η4 * vf) / (1 - η4 * vf)


    nu23_UD = E22_UD / (2G23_UD) - 1

    nu31_UD = (nu12_UD * E33_UD) / E11_UD


    peff = OrthotropicElasticParameters(;E1 = E11_UD,
                                         E2 = E22_UD, 
                                         E3 = E33_UD, 
                                         G12 = G12_UD, 
                                         G23 = G23_UD, 
                                         G31 = G31_UD, 
                                         nu12 = nu12_UD, 
                                         nu31 = nu31_UD, 
                                         nu23 = nu23_UD)

                                         
    return stiffness_matrix_voigt(peff; mandel)

end





# Adapted for Transverse Isotropic Fibers
# Ef1: Fiber Axial Modulus
# Ef2: Fiber Transverse Modulus
# Gf12: Fiber Longitudinal-Transverse Shear Modulus
# Gf23: Fiber Transverse-Transverse Shear Modulus
# nu_f12: Fiber Major Poisson's Ratio
function halpin_tsai_transverse_isotropic(Ef1, Ef2, Gf12, Gf23, nu_f12, Em, nu_m, vf, ar; mandel = true)

    # 1. Longitudinal Modulus (E11) 
    # Usually Rule of Mixtures is more accurate for UD, but Halpin-Tsai with 2*ar works too
    ζ1 = 2 * ar
    η1 = (Ef1 / Em - 1) / (Ef1 / Em + ζ1)
    E11_UD = Em * (1 + ζ1 * η1 * vf) / (1 - η1 * vf)

    # 2. Major Poisson's Ratio (nu12)
    nu12_UD = nu_f12 * vf + nu_m * (1 - vf)

    # 3. Transverse Modulus (E22 & E33)
    # CRITICAL: Use Ef2 here, not Ef1
    ζ2 = 2
    η2 = (Ef2 / Em - 1) / (Ef2 / Em + ζ2)
    E22_UD = E33_UD = Em * (1 + ζ2 * η2 * vf) / (1 - η2 * vf)

    # 4. Longitudinal Shear Modulus (G12 & G13)
    # CRITICAL: Use Gf12 here
    Gm = Em / (2 * (1 + nu_m))
    ζ3 = 1
    η3 = (Gf12 / Gm - 1) / (Gf12 / Gm + ζ3)
    G12_UD = G31_UD = Gm * (1 + ζ3 * η3 * vf) / (1 - η3 * vf)

    # 5. Transverse Shear Modulus (G23)
    ζ4 = (1 + nu_m) / (3 - nu_m - 4 * nu_m^2)
    η4 = (Gf23 / Gm - 1) / (Gf23 / Gm + ζ4)
    G23_UD = Gm * (1 + ζ4 * η4 * vf) / (1 - η4 * vf)

    # 6. Secondary Poisson's Ratios
    # Using the relationship for transverse isotropy
    nu23_UD = E22_UD / (2 * G23_UD) - 1
    nu21_UD = (nu12_UD * E22_UD) / E11_UD # Reciprocity
    nu31_UD = nu21_UD

    peff = OrthotropicElasticParameters(;E1 = E11_UD,
                                         E2 = E22_UD, 
                                         E3 = E33_UD, 
                                         G12 = G12_UD, 
                                         G23 = G23_UD, 
                                         G31 = G31_UD, 
                                         nu12 = nu12_UD, 
                                         nu31 = nu31_UD, 
                                         nu23 = nu23_UD)


    return stiffness_matrix_voigt(peff; mandel)
end


function halpin_tsai(pm::IsotropicElasticParameters, pf::TransverseIsotropicElasticParameters, vf, ar)
    Em, nu_m = pm.E_modulus, pm.nu
    Ef1, Ef2, Gf12, Gf23, nu21_f = pf.E1, pf.E2, pf.G12, pf.G23, pf.nu21 
    nu12_f = nu21_f * Ef1 / Ef2

    return halpin_tsai_transverse_isotropic(Ef1, Ef2, Gf12, Gf23, nu12_f, Em, nu_m, vf, ar)
end



#######    ORIENTATION AVERAGING        ###############



# Advani-Tucker Orientation Averaging
function orientation_average(C_aligned::AbstractMatrix, a::AbstractOrientationTensor; 
                            closure_type = HybridClosure::Type{<:AbstractClosure}, 
                            mandel = true,
                            )
    

    tens = convert_66_to_3333(C_aligned; mandel)
    Cavg = orientation_average(tens,a; closure_type)
    
    return mandel ? tomandel(Cavg) : tovoigt(Cavg)
end


function orientation_averaging_coefficients(C::AbstractMatrix; mandel = true)
    tens = convert_66_to_3333(C; mandel)
    
    return orientation_averaging_coefficients(tens)
end

function orientation_averaging_coefficients(tens::SymmetricTensor{4,dim}) where {dim}
    

    B1 = tens[1,1,1,1] + tens[2,2,2,2] - 2 * tens[1,1,2,2] - 4 * tens[1,2,1,2]
           
    B2 = tens[1,1,2,2] - tens[2,2,3,3]

    B3 = tens[1,2,1,2] + 1/2 * (tens[2,2,3,3] - tens[2,2,2,2])

    B4 = tens[2,2,3,3]

    B5 = 1/2 * (tens[2,2,2,2] - tens[2,2,3,3])

    return (B1, B2, B3, B4, B5)
end


function orientation_average(tens::SymmetricTensor{4, 3}, a::AbstractOrientationTensor; 
                            closure_type = HybridClosure::Type{<:AbstractClosure})

    
    (B1, B2, B3, B4, B5) = orientation_averaging_coefficients(tens)

    a4 = closure(a, closure_type)
    a2 = to_matrix(a)

    Cavg = SymmetricTensor{4, 3}(
        (i, j, k, l) -> 
                    B1 * a4[i, j, k, l] + 
                    B2 * (a2[i, j] * δ(k,l) + a2[k, l] * δ(i, j)) +
                    B3 * (a2[i, k] * δ(j, l) + a2[i,l] * δ(j, k) + a2[j, l] * δ(i, k) + a2[j, k] * δ(i, l)) +
                    B4 * (δ(i, j) * δ(k, l)) + 
                    B5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
                    )

    return Cavg

end