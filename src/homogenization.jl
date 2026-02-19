
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
                                        )

    Cm = stiffness_matrix_voigt(matrix_props)
    Cf = stiffness_matrix_voigt(fiber_props)
    
    if !isnothing(mass_fraction) #overwrite volume fraction
        volume_fraction = calc_vol_fraction(mass_fraction, density_fiber, density_matrix)
    end
    
    C_aligned = mori_tanaka(Cm, Cf, volume_fraction, aspect_ratio, matrix_props.nu) 
    C_averaged = orientation_average(C_aligned, orientation_tensor)

    return extract_orthotropic_constants(C_averaged)
end

# Main wrapper function to be called by users
function compute_orthotropic_properties(Em, num, Ef, nuf, vf, AR, a11, a22)

    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    C_aligned = mori_tanaka(Cm, Cf, vf, AR, num)
    C_avg = orientation_average(C_aligned, OrientationTensor(a11, a22))
    return extract_orthotropic_constants(C_avg)
end




# function eshelby_tensor_spheroid(nu_m, AR)
#     S = @MMatrix zeros(6, 6)
    
#     if AR ≈ 1.0
#         # --- SPHERICAL CASE (Isotropic) ---
#         s_diag = (7 - 5*nu_m) / (15 * (1 - nu_m))
#         s_off  = (5*nu_m - 1) / (15 * (1 - nu_m))
#         s_shear = (4 - 5*nu_m) / (15 * (1 - nu_m))
        
#         S[1,1]=S[2,2]=S[3,3] = s_diag
#         S[1,2]=S[1,3]=S[2,1]=S[2,3]=S[3,1]=S[3,2] = s_off
#         S[4,4]=S[5,5]=S[6,6] = 2 * s_shear
        
#     else
#         # --- SPHEROIDAL CASES (1-axis is the symmetry axis) ---
#         θ = AR
#         # Calculate g-factor based on shape
#         if θ > 1.0
#             # Prolate (Fibers)
#             g = (θ / (θ^2 - 1.0)^1.5) * (θ * sqrt(θ^2 - 1.0) - acosh(θ))
#         else
#             # Oblate (Disks)
#             g = (θ / (1.0 - θ^2)^1.5) * (acos(θ) - θ * sqrt(1.0 - θ^2))
#         end
        
#         # Longitudinal Component (along 1-axis)
#         S1111 = (1.0 / (2.0*(1.0-nu_m))) * (1.0 - 2.0*nu_m + (3.0*θ^2 - 1.0)/(θ^2 - 1.0) - (1.0 - 2.0*nu_m + 3.0*θ^2/(θ^2 - 1.0))*g)
        
#         # Transverse Component (along 2 and 3 axes)
#         S2222 = (3.0 / (8.0*(1.0-nu_m))) * (θ^2/(θ^2 - 1.0)) + (1.0 / (4.0*(1.0-nu_m))) * (1.0 - 2.0*nu_m - 9.0/(4.0*(θ^2 - 1.0)))*g
#         S3333 = S2222
        
#         # Coupling: Longitudinal strain due to Transverse stress
#         S1122 = (1.0 / (2.0*(1.0-nu_m))) * (-(θ^2/(θ^2 - 1.0)) + (1.0 - 2.0*nu_m + 3.0/(2.0*(θ^2 - 1.0)))*g)
#         S1133 = S1122
        
#         # Coupling: Transverse strain due to Longitudinal stress
#         S2211 = (1.0 / (4.0*(1.0-nu_m))) * (-(θ^2/(θ^2 - 1.0)) + (3.0*θ^2/(θ^2 - 1.0) - (1.0 - 2.0*nu_m))*g)
#         S3311 = S2211
        
#         # Coupling: Transverse-Transverse (2-3 coupling)
#         S2233 = (1.0 / (8.0*(1.0-nu_m))) * (θ^2/(θ^2 - 1.0) - (1.0 - 2.0*nu_m + 3.0/(4.0*(θ^2 - 1.0)))*g)
#         S3322 = S2233
        
#         # Shear Terms
#         # S2323 (Transverse plane shear)
#         S2323 = (1.0 / (8.0*(1.0-nu_m))) * (θ^2/(θ^2 - 1.0) + (1.0 - 2.0*nu_m - 3.0/(4.0*(θ^2 - 1.0)))*g)
#         # S1212 (Longitudinal-Transverse shear)
#         S1212 = (1.0 / (4.0*(1.0-nu_m))) * (1.0 - 2.0*nu_m - (θ^2 + 1.0)/(θ^2 - 1.0) - 0.5*(1.0 - 2.0*nu_m - 3.0*(θ^2 + 1.0)/(θ^2 - 1.0))*g)
#         S1313 = S1212
        
#         # Map to 6x6 Voigt Matrix (1=xx, 2=yy, 3=zz, 4=yz, 5=zx, 6=xy)
#         S[1,1] = S1111
#         S[2,2] = S2222
#         S[3,3] = S3333
        
#         S[1,2] = S[1,3] = S1122
#         S[2,1] = S[3,1] = S2211
#         S[2,3] = S[3,2] = S2233
        
#         # Multiply by 2 for engineering shear strain mapping in Voigt form
#         S[4,4] = 2.0 * S2323
#         S[5,5] = 2.0 * S1313
#         S[6,6] = 2.0 * S1212
#     end
    
#     return SMatrix{6,6}(S...)
# end

"""
    eshelby_spheroid(ar, nu)
Calculates the Eshelby tensor S for a prolate spheroid (aspect ratio ar > 1)
aligned along the 1-axis in an isotropic matrix with Poisson's ratio nu.
"""
function eshelby_tensor_spheroid(nu, ar)
    # Auxiliary variable g for prolate spheroids
    # ar = length / diameter
    g = (ar / (ar^2 - 1)^1.5) * (ar * sqrt(ar^2 - 1) - acosh(ar))
    
    # Components S_ijkl in Voigt notation
    S = MMatrix{6,6}(zeros(6, 6))
    
    # Denominator factor
    f = 8 * π * (1 - nu)
    
    # S1111 (S[1,1])
    S[1,1] = (1 / f) * (4π * (1 - 2nu) + 3 * 4π * (ar^2 / (ar^2 - 1)) * (1 - g))
    
    # S2222 = S3333 (S[2,2], S[3,3])
    S[2,2] = S[3,3] = (3 / (8 * (1 - nu))) * (ar^2 / (ar^2 - 1)^2) + (1 / (4 * (1 - nu))) * (1 - 2nu + (9 / (4 * (ar^2 - 1)))) * g
    # Simplified standard forms for S-tensor components:
    I1 = 4π * (1 - g) * (ar^2 / (ar^2 - 1))
    I2 = 2π * g * (ar / (ar^2 - 1)^0.5) # This is a simplified branch
    
    # Using Tandon & Weng (1984) analytical expressions for aligned fibers:
    # We'll use the precise components for a transverse isotropic S-tensor
    S11 = (1 / (2 * (1 - nu))) * (1 - 2nu + (3ar^2 - 1)/(ar^2 - 1) - (1 - 2nu + (3ar^2)/(ar^2 - 1)) * g)
    S22 = S33 = (3 / (8 * (1 - nu))) * (ar^2 / (ar^2 - 1)) - (1 / (4 * (1 - nu))) * (1 - 2nu + (9 / (4 * (ar^2 - 1)))) * g
    S23 = S32 = (1 / (4 * (1 - nu))) * (ar^2 / (2 * (ar^2 - 1)) - (1 - 2nu + (3 / (4 * (ar^2 - 1)))) * g)
    S21 = S31 = -(1 / (2 * (1 - nu))) * (ar^2 / (ar^2 - 1)) + (1 / (4 * (1 - nu))) * (3ar^2 / (ar^2 - 1) - (1 - 2nu)) * g
    S12 = S13 = -(1 / (2 * (1 - nu))) * (1 - 2nu + 1/(ar^2 - 1)) + (1 / (2 * (1 - nu))) * (1 - 2nu + 3/(2 * (ar^2 - 1))) * g
    S44 = (S22 - S23) / 2
    S55 = S66 = (1 / (4 * (1 - nu))) * (1 - 2nu - (ar^2 + 1)/(ar^2 - 1) + (1/2 * (1 - 2nu + (3 * (ar^2 + 1))/(ar^2 - 1))) * g)

    S[1,1] = S11
    S[2,2] = S[3,3] = S22
    S[2,3] = S[3,2] = S23
    S[2,1] = S[3,1] = S21
    S[1,2] = S[1,3] = S12
    S[4,4] = S44
    S[5,5] = S[6,6] = S55
    
    return SMatrix{6,6}(S)
end


"""
    mori_tanaka_homogenization(Cm, Ci, S, f)
Cm: Matrix Stiffness, Ci: Inclusion Stiffness, S: Eshelby, f: Volume fraction
"""
function mori_tanaka(Cm::AbstractMatrix, Cf::AbstractMatrix, vf, AR, nu_m)
    T = eltype(Cm)
    I = SMatrix{6,6, T}(LinearAlgebra.I)
    S = eshelby_tensor_spheroid(nu_m, AR)
    
    # Dilute concentration tensor (T)
    T = inv(I + S * (inv(Cm) * (Cf - Cm)))
    
    # Mori-Tanaka concentration tensor (A)
    A = T * inv((1 - vf) * I + vf * T)
    
    # Effective Stiffness
    C_eff = Cm + vf * (Cf - Cm) * A
    return SMatrix{6,6}(C_eff)
end


# function mori_tanaka(Cm::AbstractMatrix, Cf::AbstractMatrix, vf, AR, nu_m)
#     T = eltype(Cm)
#     I = SMatrix{6,6, T}(LinearAlgebra.I)
    
#     S = eshelby_tensor_spheroid(nu_m, AR)
#     # 1. Calculate the Dilute Concentration Tensor (A_dilute)
#     # A_dilute = [I + S * inv(Cm) * (Cf - Cm)]^-1
#     # This accounts for the strain in a single fiber relative to the matrix
    
#     # Note: inv(Cm) * (Cf - Cm) is effectively the difference in compliance
#     A_dilute = inv(I + S * (inv(Cm) * (Cf - Cm)))
    
#     # 2. Apply the Mori-Tanaka interaction term
#     # This accounts for the "crowding" of fibers as f increases
#     term1 = vf * (Cf - Cm) * A_dilute
#     term2 = inv((1 - vf) * I + vf * A_dilute)
    
#     C_eff = Cm + term1 * term2
    
#     return SMatrix{6,6}(C_eff)
# end

#kron delta
@inline function δ(i, j) 
    i == j ? 1 : 0
end



function I4(T=Float64)

    I4 = zeros(T, 3,3,3,3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        I4[i,j,k,l] = 1/2 *(δ(i, j) * δ(j, l) + δ(i, l) * δ(j, k))
    end
end


function convert_voigt_to_tensor(C66)

    lookup = (
        (1, 6, 5),
        (6, 2, 4),
        (5, 4, 3)
    )

    # Generate the 4th order tensor using a comprehension
    # This creates an SArray{Tuple{3,3,3,3}}
    C4 = @SArray [C66[lookup[i][j], lookup[k][l]] for i=1:3, j=1:3, k=1:3, l=1:3]
    
    return C4
end


"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
# function hybrid_closure(orientation_tensor)
#     a = orientation_tensor

#     # 1. Calculate the interaction scalar 'f'
#     # f = 1 - 27*det(a)
#     # For a 3x3 orientation tensor, det(a) ranges from 0 (aligned) to 1/27 (random)
#     f = 1.0 - 27.0 * det(a)
    
#     # Ensure f stays within physical bounds due to numerical precision
#     f = clamp(f, 0.0, 1.0)
    
#     A_lin    = @MArray zeros(3, 3, 3, 3)
#     A_quad   = @MArray zeros(3, 3, 3, 3)
#     A_hybrid = @MArray zeros(3, 3, 3, 3)
    
#     δ(i, j) = (i == j) ? 1.0 : 0.0

#     for i in 1:3, j in 1:3, k in 1:3, l in 1:3
#         # Linear Component
#         term1 = (δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) / 24.0

#         term2 = (a[i,j]*δ(k,l) + a[i,k]*δ(j,l) + a[i,l]*δ(j,k) + 
#                  a[k,l]*δ(i,j) + a[j,l]*δ(i,k) + a[j,k]*δ(i,l)) / 6.0
        
#         A_lin[i,j,k,l] = -term1 + term2
        
#         # Quadratic Component
#         A_quad[i,j,k,l] = a[i,j] * a[k,l]
        
#         # Hybrid Interpolation
#         A_hybrid[i,j,k,l] = (1.0 - f) * A_lin[i,j,k,l] + f * A_quad[i,j,k,l]
#     end
    
#     return SArray{Tuple{3,3,3,3}}(A_hybrid)
# end


function linear_closure(a)
    T = eltype(a)
    a4 = MArray{Tuple{3,3,3,3}}(zeros(T, 3, 3, 3, 3))

    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        # Linear Component
        term1 = (δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) / 24.0

        term2 = (a[i,j]*δ(k,l) + a[i,k]*δ(j,l) + a[i,l]*δ(j,k) + 
                 a[k,l]*δ(i,j) + a[j,l]*δ(i,k) + a[j,k]*δ(i,l)) / 6.0
        
        a4[i,j,k,l] = -term1 + term2
    end
    return SArray{Tuple{3,3,3,3}}(a4)
end

function quadratic_closure(a2)
    T = eltype(a2)
    a4 = MArray{Tuple{3,3,3,3}}(zeros(T, 3, 3, 3, 3))

    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        a4[i, j, k, l] = a2[i, j] * a2[k, l]
    end
    return SArray{Tuple{3,3,3,3}}(a4)
end


function hybrid_closure(a2)
    f = -0.5
    for i in 1:3
        for j in 1:3
            f += 1.5 * a2[i, j] * a2[j, i]
        end
    end
    a4_l = linear_closure(a2)
    a4_q = quadratic_closure(a2)

    return (1 - f) * a4_l + f * a4_q
end


# Advani-Tucker Orientation Averaging
function orientation_average(C_aligned, a::OrientationTensor)
    a_mat = to_matrix(a)
    A4 = hybrid_closure(a_mat)
    
    # Extract Invariants from C_aligned (Assumes Axis 1 is longitudinal)
    B1 = C_aligned[1,1] + C_aligned[3,3] - 2*C_aligned[2,3] - 4*C_aligned[5,5]
    B2 = C_aligned[2,3] - C_aligned[1,2]
    B3 = C_aligned[5,5] + 0.5*(C_aligned[1, 1] - C_aligned[1,2])
    B4 = C_aligned[1, 2]
    B5 = 0.5*(C_aligned[1,1] - C_aligned[1, 2])
    
    C_avg = @MMatrix zeros(6,6)
    # Voigt Mapping: 1->(1,1), 2->(2,2), 3->(3,3), 4->(2,3), 5->(1,3), 6->(1,2)
    v = [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]
    
    # δ(i, j) = (i == j) ? 1.0 : 0.0

    for r in 1:6, c in 1:6
        i,j = v[r]; k,l = v[c]
        
        # This is the general expression for <C>ijkl
        val = B1 * A4[i,j,k,l] + 
              B2 * (a_mat[i,j]*δ(k,l) + a_mat[k,l]*δ(i,j)) + 
              B3 * (a_mat[i,k]*δ(j,l) + a_mat[i,l]*δ(j,k) +  a_mat[j,k]*δ(i,l)) + 
              B4 * (δ(i,j)*δ(k,l)) + 
              B5 * (δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
        
        C_avg[r,c] = val
    end
    return SMatrix{6,6}(C_avg)
end



