
"""
    compute_orthotropic_properties(matrix_props::IsotropicElasticParameters,
                                        fiber_props::AbstractElasticParameters;
                                        volume_fraction = 0.2,
                                        aspect_ratio = 10.0,
                                        orientation_tensor=OrientationTensor(0.7, 0.3),
                                        )

Effective stiffness matrix (Voigt 6x6 form) for matrix and fiber
"""
function compute_orthotropic_properties(matrix_props::IsotropicElasticParameters,
                                        fiber_props::AbstractElasticParameters;
                                        volume_fraction = 0.2,
                                        mass_fraction = nothing,
                                        density_fiber = 1.0,
                                        density_matrix = 1.0,
                                        aspect_ratio = 10.0,
                                        orientation_tensor=OrientationTensor(0.7, 0.3),
                                        shape = SpheroidalInclusion,
                                        mandel = false,
                                        )

    Cm = stiffness_matrix_voigt(matrix_props)
    Cf = stiffness_matrix_voigt(fiber_props)
    
    if !isnothing(mass_fraction) #overwrite volume fraction
        volume_fraction = calc_vol_fraction(mass_fraction, density_fiber, density_matrix)
    end
    
    C_aligned = mori_tanaka(Cm, Cf, volume_fraction, aspect_ratio, matrix_props.nu; mandel) 
    C_averaged = orientation_average(C_aligned, orientation_tensor; mandel)

    # return extract_orthotropic_constants(C_averaged)
    return C_averaged
end

# Main wrapper function to be called by users
function compute_orthotropic_properties(Em, num, Ef, nuf, vf, AR, a11, a22)

    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    C_aligned = mori_tanaka(Cm, Cf, vf, AR, num)
    a = OrientationTensor(a11, a22)
    
    C_avg = orientation_average(C_aligned, a)
    # return extract_orthotropic_constants(C_avg)
end




function mt_dilute_tensor(Cm, Cf, S)
    
    I = SMatrix{6,6}(LinearAlgebra.I)
    # 1. Calculate the Dilute Concentration Tensor (A_dilute)
    # A_dilute = [I + S * inv(Cm) * (Cf - Cm)]^-1
    # This accounts for the strain in a single fiber relative to the matrix
    return inv(I + S * (inv(Cm) * (Cf - Cm)))
end

function mori_tanaka(Cm::AbstractMatrix, Cf::AbstractMatrix, vf, AR, nu_m; 
                    fiber_shape = SpheroidalInclusion, 
                    mandel = false,
                    )

    T = eltype(Cm)
    I = SMatrix{6,6, T}(LinearAlgebra.I)
    
    S= convert_3333_to_66(eshelby_tensor(fiber_shape(nu_m, AR)); mandel)
    
    A_dilute = mt_dilute_tensor(Cm, Cf, S)
    
    # 2. Apply the Mori-Tanaka interaction term
    # This accounts for the "crowding" of fibers as f increases
    term1 = vf * (Cf - Cm) * A_dilute
    term2 = inv((1 - vf) * I + vf * A_dilute)
    
    C_eff = Cm + term1 * term2
    
    return SMatrix{6,6}(C_eff)
end

function mori_tanaka(pm::IsotropicElasticParameters, fibers; mandel = false)

    Cm = stiffness_matrix_voigt(pm; mandel)
    invCm = inv(Cm)

    T  = eltype(Cm)
    MT = SMatrix{6,6, T}
    
    nu_m = pm.nu


    #precompute everything
    Cfs =      Vector{MT}()
    Adils =    Vector{MT}()
    # Stensors = Vector{MT}()
    vfs = (fiber.volume_fraction for fiber in fibers)
    for fiber in fibers
        vf = fiber.volume_fraction
        Cf = stiffness_matrix_voigt(fiber.elastic_properties; mandel)
        push!(Cfs, Cf)
        S = eshelby_tensor(fiber.shape(nu_m, fiber.AR))
        # push!(Stensors, S)
        Adil = inv(vf * I + S * invCm * (Cf - Cm))
        push!(Adils, Adil)
    end

    Ceff = MMatrix{6,6}(Cm)
    for (i, fiber) in enumerate(fibers)
        
        Adil = Adils[i]
        vf = fiber.volume_fraction
        Cf = Cfs[i]
        ArMT = Adil * inv(vf * I + sum(vfs[i] * Adils[i] for i in eachindex(fibers)))
        Ceff += vf * (Cf - Cm) * ArMT
    end
    return SMatrix{6,6}(C_eff)
end





function halpin_tsai(Ef, Em, nu_f, nu_m, vf, ar)

    ζ1 = 2ar
    η1 = (Ef/Em - 1)/ (Ef/Em + ζ1)

    E11_UD = Em * (1 + ζ1 * η1 * vf) / (1 - η1 * vf)
    @info "E11 $E11_UD"
    nu12_UD = nu_f * vf + nu_m * (1 - vf)

    ζ2 = 2
    η2 = (Ef/Em - 1) / (Ef/Em + ζ2)

    E22_UD = E33_UD = Em * (1 + ζ2 * η2 * vf) / (1 - η2 * vf)
    @info "E22 $E22_UD"
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

    S11_UD = 1 / E11_UD
    S22_UD = 1 / E22_UD
    S33_UD = 1 / E33_UD
    S12_UD = -nu12_UD / E11_UD
    S13_UD = -nu31_UD / E33_UD
    S23_UD = -nu23_UD / E22_UD
    S44_UD = 1/G23_UD
    S55_UD = 1/G31_UD
    S66_UD = 1/G12_UD

    S_UD = @SMatrix [S11_UD  S12_UD  S13_UD  0      0     0;
                     S12_UD  S22_UD  S23_UD  0      0     0;
                     S13_UD  S23_UD  S33_UD  0      0     0;
                       0       0       0    S44_UD  0     0;
                       0       0       0     0   S55_UD   0;
                       0       0       0     0      0  S66_UD]  
    return (inv(S_UD), (;E11_UD, E22_UD, G12_UD, G23_UD, nu12_UD, nu23_UD))

end






# Advani-Tucker Orientation Averaging
function orientation_average(C_aligned, a::AbstractOrientationTensor; 
                            closure_type = HybridClosure::Type{<:AbstractClosure}, 
                            mandel = false,
                            )
    
    (B1, B2, B3, B4, B5) = orientation_averaging_coefficients(C_aligned) #this works
    a2 = to_matrix(a)
    a4 = closure(a, closure_type)

    # The Advani-Tucker orientation averaging formula:
#     # <C_ijkl> = b1(a_ijkl) + b2(a_ij*d_kl + a_kl*d_ij) + 
#     #            b3(a_ik*d_jl + a_il*d_jk + a_jl*d_ik + a_jk*d_il) + 
#     #            b4(d_ij*d_kl) + b5(d_ik*d_jl + d_il*d_jk)

    Cavg = SArray{Tuple{3,3,3,3}}( 
                    B1 * a4[i, j, k, l] + 
                    B2 * (a2[i, j] * δ(k,l) + a2[k, l] * δ(i, j)) +
                    B3 * (a2[i, k] * δ(j, l) + a2[i,l] * δ(j, k) + a2[j, l] * δ(i, k) + a2[j, k] * δ(i, l)) +
                    B4 * (δ(i, j) * δ(k, l)) + 
                    B5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
                    for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

    return convert_3333_to_66(Cavg; mandel)
end

function orientation_averaging_coefficients(C)
    # Define the 5 invariants for a transversely isotropic material 
    # based on the aligned stiffness tensor components.
    # These are derived from the terms in Advani & Tucker (1987).
    # Extract components from the aligned tensor (assuming 1 is the fiber direction)
#     b1 = Caligned[1,1,1,1] + Caligned[2,2,2,2] - 2*Caligned[1,1,2,2] - 4*Caligned[1,2,1,2]
#     b2 = Caligned[1,1,2,2] - Caligned[2,2,3,3]
#     b3 = Caligned[1,2,1,2] + 0.5 * (Caligned[2,2,3,3] - Caligned[2,2,2,2])
#     b4 = Caligned[2,2,3,3]
#     b5 = 0.5 * (Caligned[2,2,2,2] - Caligned[2,2,3,3])

    B1 = C[1,1] + C[2,2]  - 2C[1, 2] - 4C[6,6]
    B2 = C[1,2] - C[2,3]
    B3 = C[6,6] + 1/2 * (C[2,3] - C[2,2])
    B4 = C[2,3]
    B5 = 1/2 * (C[2,2] - C[2,3])

    return (B1, B2, B3, B4, B5)
end

using LinearAlgebra


#Gemini generated
# """
#     orientation_average(Caligned, N2, N4)

# Computes the orientation-averaged stiffness tensor for SFRP.
# Assumes Caligned is a 3x3x3x3 tensor (or a 6x6 Voigt matrix converted to 4D).
# N2 is the 2nd order orientation tensor (3x3).
# N4 is the 4th order orientation tensor (3x3x3x3).
# """
# function orientation_average(Caligned, N2, N4)
#     # Define the 5 invariants for a transversely isotropic material 
#     # based on the aligned stiffness tensor components.
#     # These are derived from the terms in Advani & Tucker (1987).
    
#     # Extract components from the aligned tensor (assuming 1 is the fiber direction)
#     b1 = Caligned[1,1,1,1] + Caligned[2,2,2,2] - 2*Caligned[1,1,2,2] - 4*Caligned[1,2,1,2]
#     b2 = Caligned[1,1,2,2] - Caligned[2,2,3,3]
#     b3 = Caligned[1,2,1,2] + 0.5 * (Caligned[2,2,3,3] - Caligned[2,2,2,2])
#     b4 = Caligned[2,2,3,3]
#     b5 = 0.5 * (Caligned[2,2,2,2] - Caligned[2,2,3,3])

#     C_avg = zeros(3, 3, 3, 3)

#     # The Advani-Tucker orientation averaging formula:
#     # <C_ijkl> = b1(a_ijkl) + b2(a_ij*d_kl + a_kl*d_ij) + 
#     #            b3(a_ik*d_jl + a_il*d_jk + a_jl*d_ik + a_jk*d_il) + 
#     #            b4(d_ij*d_kl) + b5(d_ik*d_jl + d_il*d_jk)

#     δ = I(3) # Kronecker Delta

#     for i in 1:3, j in 1:3, k in 1:3, l in 1:3
#         # Term 1: 4th order orientation
#         term1 = b1 * N4[i,j,k,l]
        
#         # Term 2: Coupling 2nd order with identity
#         term2 = b2 * (N2[i,j]*δ[k,l] + N2[k,l]*δ[i,j])
        
#         # Term 3: Symmetric coupling
#         term3 = b3 * (N2[i,k]*δ[j,l] + N2[i,l]*δ[j,k] + N2[j,l]*δ[i,k] + N2[j,k]*δ[i,l])
        
#         # Term 4 & 5: Isotropic components
#         term4 = b4 * (δ[i,j]*δ[k,l])
#         term5 = b5 * (δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k])
        
#         C_avg[i,j,k,l] = term1 + term2 + term3 + term4 + term5
#     end

#     return C_avg
# end