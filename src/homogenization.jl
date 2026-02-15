
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
                                        aspect_ratio = 10.0,
                                        orientation_tensor=OrientationTensor(0.7, 0.3),
                                        )

    Cm = stiffness_matrix_voigt(matrix_props)
    Cf = stiffness_matrix_voigt(fiber_props)
    
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


# Full Analytical Eshelby Tensor for Prolate Spheroids
function eshelby_tensor_prolate(nu_m, AR)
    S = @MMatrix zeros(6, 6)
    # Numerical safety for spheres
    θ = AR <= 1.0 ? 1.0001 : AR
    g = (θ / (θ^2 - 1)^1.5) * (θ * sqrt(θ^2 - 1) - acosh(θ))
    
    S1111 = (1 / (2*(1-nu_m))) * (1 - 2nu_m + (3θ^2 - 1)/(θ^2 - 1) - (1 - 2nu_m + 3θ^2/(θ^2-1))*g)
    S2222 = (3 / (8*(1-nu_m))) * (θ^2/(θ^2-1)) + (1 / (4*(1-nu_m))) * (1 - 2nu_m - 9/(4*(θ^2-1)))*g
    S3333 = S2222
    S1122 = (1 / (2*(1-nu_m))) * (- (θ^2/(θ^2-1)) + (1 - 2nu_m + 3/(2*(θ^2-1)))*g)
    S2233 = (1 / (8*(1-nu_m))) * (θ^2/(θ^2-1) - (1 - 2nu_m + 3/(4*(θ^2-1)))*g)
    S2211 = (1 / (4*(1-nu_m))) * (- (θ^2/(θ^2-1)) + (3θ^2/(θ^2-1) - (1-2nu_m))*g)
    S1212 = (1 / (4*(1-nu_m))) * (1 - 2nu_m - (θ^2+1)/(θ^2-1) - 0.5*(1 - 2nu_m - 3*(θ^2+1)/(θ^2-1))*g)
    S2323 = (1 / (8*(1-nu_m))) * (θ^2/(θ^2-1) + (1 - 2nu_m - 3/(4*(θ^2-1)))*g)
    
    S[1,1]=S1111; S[2,2]=S3333; S[3,3]=S3333
    S[1,2]=S[1,3]=S1122; S[2,1]=S[3,1]=S2211; S[2,3]=S[3,2]=S2233
    S[4,4]=2S2323; S[5,5]=S[6,6]=2S1212
    return SMatrix{6,6}(S...)
end

# Mori-Tanaka Aligned Homogenization
function mori_tanaka(Cm, Cf, vf, AR, nu_m)
    I6 = Matrix{eltype(Cm)}(I, 6, 6)
    S = eshelby_tensor_prolate(nu_m, AR)
    #dilute concentration tensor
    A_dil = inv(I6 + S * (inv(Cm) * (Cf - Cm)))
    A_MT = A_dil * inv((1 - vf) * I6 + vf * A_dil)
    return Cm + vf * (Cf - Cm) * A_MT
end


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
function hybrid_closure(orientation_tensor)
    a = orientation_tensor

    # 1. Calculate the interaction scalar 'f'
    # f = 1 - 27*det(a)
    # For a 3x3 orientation tensor, det(a) ranges from 0 (aligned) to 1/27 (random)
    f = 1.0 - 27.0 * det(a)
    
    # Ensure f stays within physical bounds due to numerical precision
    f = clamp(f, 0.0, 1.0)
    
    A_lin    = @MArray zeros(3, 3, 3, 3)
    A_quad   = @MArray zeros(3, 3, 3, 3)
    A_hybrid = @MArray zeros(3, 3, 3, 3)
    
    δ(i, j) = (i == j) ? 1.0 : 0.0

    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        # Linear Component
        term1 = (δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) / 24.0

        term2 = (a[i,j]*δ(k,l) + a[i,k]*δ(j,l) + a[i,l]*δ(j,k) + 
                 a[k,l]*δ(i,j) + a[j,l]*δ(i,k) + a[j,k]*δ(i,l)) / 6.0
        
        A_lin[i,j,k,l] = -term1 + term2
        
        # Quadratic Component
        A_quad[i,j,k,l] = a[i,j] * a[k,l]
        
        # Hybrid Interpolation
        A_hybrid[i,j,k,l] = (1.0 - f) * A_lin[i,j,k,l] + f * A_quad[i,j,k,l]
    end
    
    return SArray{Tuple{3,3,3,3}}(A_hybrid)
end


# Advani-Tucker Orientation Averaging
function orientation_average(C_aligned, a::OrientationTensor)
    a11, a22 = a.a11, a.a22
    a_mat = to_matrix(OrientationTensor(a11, a22))
    A4 = hybrid_closure(a_mat)
    
    # Transversely Isotropic Invariants
    B1 = C_aligned[1,1] + C_aligned[2,2] - 2*C_aligned[1,2] - 4*C_aligned[6,6]
    B2 = C_aligned[1,2] - C_aligned[2,3]
    B3 = C_aligned[6,6] + 0.5*(C_aligned[2,3] - C_aligned[2,2])
    B4 = C_aligned[2,3]
    B5 = 0.5*(C_aligned[2,2] - C_aligned[2,3])
    
    C_avg = @MMatrix zeros(6,6)
    v = @SVector [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]

    
    for r in 1:6, c in 1:6
        i,j = v[r]; k,l = v[c]
        C_avg[r,c] = B1*A4[i,j,k,l] + B2*(a_mat[i,j]*δ(k,l) + a_mat[k,l]*δ(i,j)) + 
                     B3*(a_mat[i,k]*δ(j,l) + a_mat[i,l]*δ(j,k) + a_mat[j,l]*δ(i,k) + a_mat[j,k]*δ(i,l)) + 
                     B4*(δ(i,j)*δ(k,l)) + B5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    end
    return SMatrix{6,6}(C_avg)
end




