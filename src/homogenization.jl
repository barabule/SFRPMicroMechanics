
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
    N2 = to_matrix(OrientationTensor(a11, a22))
    C_avg = orientation_average(C_aligned, N2)
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
    
    # S = convert_3333_to_66(eshelby_tensor(fiber_shape(nu_m, AR)); mandel = false)
    S= convert_3333_to_66(eshelby_tensor(fiber_shape(nu_m, AR)); mandel)
    
    A_dilute = mt_dilute_tensor(Cm, Cf, S)
    
    # 2. Apply the Mori-Tanaka interaction term
    # This accounts for the "crowding" of fibers as f increases
    term1 = vf * (Cf - Cm) * A_dilute
    term2 = inv((1 - vf) * I + vf * A_dilute)
    
    C_eff = Cm + term1 * term2
    
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
function orientation_average(C_aligned, a2; closure = hybrid_closure, mandel = false)
    
    (B1, B2, B3, B4, B5) = orientation_averaging_coefficients(C_aligned) #this works
    
    a4 = closure(a2)

    Cavg = SArray{Tuple{3,3,3,3}}( B1 * a4[i, j, k, l] + 
                                   B2 * (a2[i, j] * δ(k,l) + a2[k, l] * δ(i, j)) +
                                   B3 * (a2[i, k] * δ(j, l) + a2[j, l] * δ(i, k) + a2[i,l] * δ(j, k) + a2[j, k] * δ(i, l)) +
                                   B4 * (δ(i, j) * δ(k, l)) + 
                                   B5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
                                   for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

    return convert_3333_to_66(Cavg; mandel)
end





function orientation_averaging_coefficients(C)
    
    B1 = C[1,1] + C[2,2]  - 2C[1, 2] - 4C[6,6]
    B2 = C[1,2] - C[2,3]
    B3 = C[6,6] + 1/2 * (C[2,3] - C[2,2])
    B4 = C[2,3]
    B5 = 1/2 * (C[2,2] - C[2,3])

    return (B1, B2, B3, B4, B5)
end

