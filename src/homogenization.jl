#kron delta
@inline function δ(i, j) 
    i == j ? 1 : 0
end




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


function eshelby_tensor_spheroid(nu::T, ar::T) where T

    

    if ar ≈ 1 #sphere 
        fac1 = 15 * (1 - nu)
        
        S = SArray{Tuple{3,3,3,3}, T}( (1 / fac1 * (5nu - 1) * (δ(i, j) * δ(k, l)) + 
                    (4 - 5nu) * (δ(i,k) * δ(j,l) + δ(i,l) * δ(j,k)))
                    for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
        
    else

        if ar > 1.0 # Prolate (Fibers)
        
            g = (ar / sqrt((ar^2 - 1.0)^3)) * (ar * sqrt(ar^2 - 1.0) - acosh(ar))
        else
            # Oblate (Disks)
            g = (ar / sqrt((1.0 - ar^2)^3)) * (acos(ar) - ar * sqrt(1 - ar^2))
        end
        
        S = MArray{Tuple{3, 3, 3, 3},T}(zeros(3,3,3,3))

        S[1,1,1,1] = 1 / (2 * (1 - nu)) * ( 1 - 2nu + (3ar^2 - 1) / (ar^2 - 1) - (1 - 2nu + 3ar / (ar^2 -1)) * g)

        S[2,2,2,2] = S[3,3,3,3] = (3 / (8 * (1 - nu)) * ar^2 / (ar^2 -1)) * g 

        S[2,2,3,3] = S[3,3,2,2] = 1 / (4 * (1 - nu)) * (ar^2 / (ar^2 -1)) - (1 - 2nu + 3 / (4 * (ar^2 -1))) * g

        S[2,2,1,1] = S[3,3,1,1] = -1  / (2 * (1 - nu)) * ar^2 / (ar^2 - 1) + 1 / (4 * (1 - nu)) * (3 * ar^2 / (ar^2 - 1) - (1 - 2nu)) * g

        S[1, 1, 2, 2] = S[1, 1, 3, 3] = -1 / (2 * (1 - nu)) * (1 - 2nu + 1 / (ar^2 -1)) + 1/(2 * (1 - nu)) * (1 - 2nu + 3/(2 * (ar^2 -1))) * g

        S[2, 3, 2, 3] = S[3, 2, 3, 2] = S[3, 2, 2, 3] = S[2, 3, 3, 2] = 1 / (4 * (1 - nu)) * (ar^2 / (2 * (ar^2 - 1)))

        S[1,2,1,2 ] = S[1,3,1,3] = S[2,1,2,1] = S[2,1,1,2] = S[1,2,2,1] = S[3,1,3,1] = S[3,1,1,3] = S[1,3,3,1] =
            1 / (4 * (1 - nu)) * (1 - 2nu - (ar^2 + 1) / (ar^2 -1) - 1/2 * (1 - 2nu - 3 * (ar^2 + 1) / (ar^2 -1)) * g)

    end

    

    return tens2eng(S; mandel = true)

end




function mori_tanaka(Cm::AbstractMatrix, Cf::AbstractMatrix, vf, AR, nu_m)
    T = eltype(Cm)
    I = SMatrix{6,6, T}(LinearAlgebra.I)
    
    S = eshelby_tensor_spheroid(nu_m, AR)
    # 1. Calculate the Dilute Concentration Tensor (A_dilute)
    # A_dilute = [I + S * inv(Cm) * (Cf - Cm)]^-1
    # This accounts for the strain in a single fiber relative to the matrix
    
    # Note: inv(Cm) * (Cf - Cm) is effectively the difference in compliance
    A_dilute = inv(I + S * (inv(Cm) * (Cf - Cm)))
    
    # 2. Apply the Mori-Tanaka interaction term
    # This accounts for the "crowding" of fibers as f increases
    term1 = vf * (Cf - Cm) * A_dilute
    term2 = inv((1 - vf) * I + vf * A_dilute)
    
    C_eff = Cm + term1 * term2
    
    return SMatrix{6,6}(C_eff)
end





function I4(T=Float64)

    I4 = zeros(T, 3,3,3,3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        I4[i,j,k,l] = 1/2 *(δ(i, j) * δ(j, l) + δ(i, l) * δ(j, k))
    end
end








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

"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
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



