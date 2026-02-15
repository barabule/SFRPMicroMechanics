


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

@inline function δ(i, j) 
    i == j ? 1 : 0
end

function I4(T=Float64)

    I4 = zeros(T, 3,3,3,3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        I4[i,j,k,l] = 1/2 *(δ(i, j) * δ(j, l) + δ(i, l) * δ(j, k))
    end
end


# # This is equivalent to the math: σ_ij = C_ijkl * ε_kl
# @tullio σ[i, j] := C4[i, j, k, l] * ε[k, l]



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

function rotate_tensor(C, R)
    @tullio C_rot[i,j,k,l] := R[i,m] * R[j,n] * R[k,o] * R[l,p] * C[m,n,o,p]
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
function orientation_average(C_aligned, a11, a22)
    
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


# Main wrapper function to be called by users
function compute_orthotropic_properties(Em, num, Ef, nuf, vf, AR, a11, a22)

    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    C_aligned = mori_tanaka(Cm, Cf, vf, AR, num)
    C_avg = orientation_average(C_aligned, a11, a22)
    return extract_orthotropic_constants(C_avg)
end



"""
    compute_sfrp_cte(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22)
Returns [alpha1, alpha2, alpha3] in the principal material directions.
"""
function compute_sfrp_cte(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22)
    # 1. Setup Stiffness
    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    I6 = Matrix{Float64}(I, 6, 6)
    S_eshelby = eshelby_tensor_prolate(num, AR)
    
    # 2. MT Concentration Tensor
    Adil = inv(I6 + S_eshelby * (inv(Cm) * (Cf - Cm)))
    Amt = Adil * inv((1 - vf) * I6 + vf * Adil)
    
    # 3. Aligned Effective Stiffness
    C_aligned = Cm + vf * (Cf - Cm) * Amt
    
    # 4. Aligned CTE (Rosen and Hashin logic)
    # Transformation to Voigt vector [a1, a2, a3, 0, 0, 0]
    alpha_m_vec = [alpham, alpham, alpham, 0, 0, 0]
    alpha_f_vec = [alphaf, alphaf, alphaf, 0, 0, 0]
    
    # Aligned CTE vector
    # alpha_eff = alpha_m + vf * inv(C_aligned) * Cf * Amt * (alpha_f - alpha_m)
    term = vf * inv(C_aligned) * (Cf * Amt * (alpha_f_vec .- alpha_m_vec))
    alpha_aligned_vec = alpha_m_vec .+ term

    # 5. Orientation Averaging
    # For CTE (a 2nd order tensor), averaging is simpler than stiffness:
    # alpha_ij = (a1_aligned - a2_aligned) * a_ij + a2_aligned * delta_ij
    a1 = alpha_aligned_vec[1]
    a2 = alpha_aligned_vec[2] # Transverse aligned CTE
    
    a33 = 1.0 - a11 - a22
    alpha1 = (a1 - a2) * a11 + a2
    alpha2 = (a1 - a2) * a22 + a2
    alpha3 = (a1 - a2) * a33 + a2
    
    return [alpha1, alpha2, alpha3]
end

# """
# # 30% GF Nylon 66
# res_cte = compute_sfrp_cte(2800.0, 0.35, 80e-6,  # Matrix
# 72000.0, 0.22, 5e-6, # Fiber
# 0.18, 20.0,          # vf and AR
# 0.75, 0.15)          # Orientation

# println("CTE-1 (Flow):      $(round(res_cte[1]*1e6, digits=2)) ppm/°C")
# println("CTE-2 (Crossflow): $(round(res_cte[2]*1e6, digits=2)) ppm/°C")
# println("CTE-3 (Thickness): $(round(res_cte[3]*1e6, digits=2)) ppm/°C")
# """
"""
estimate_as_molded_ar(L_initial_mm, diameter_um, vf, processing_severity)

Estimates the aspect ratio after injection molding.
- processing_severity: 0.1 (gentle) to 1.0 (aggressive/thin gates)
- Typical L_limit for glass is ~150um.
"""
function estimate_as_molded_ar(L_initial_mm, diameter_um, vf, severity=0.5)
    L_initial_um = L_initial_mm * 1000.0
    
    # L_limit is the length below which fibers rarely break (saturation)
    # This is typically 15-20x the diameter
    L_limit = diameter_um * 15.0 
    
    # Breaking rate increases with fiber concentration (vf) and severity
    # We use an exponential decay model
    decay_constant = severity * (1.0 + 2.0 * vf)
    
    L_final = L_limit + (L_initial_um - L_limit) * exp(-decay_constant)
    
    return L_final / diameter_um
end


#### MULTI FIBER
"""
compute_hybrid_sfrp(matrix_props, fiber_list, a11, a22)
'fiber_list' is an array of dicts: [vf, E, nu, AR]
"""
function compute_hybrid_sfrp(E_m, nu_m, fibers, a11, a22)
    Cm = isotropic_stiffness(E_m, nu_m)
    I6 = Matrix{Float64}(I, 6, 6)
    
    # 1. Calculate Dilute Tensors for each fiber species
    A_dil_list = []
    v_total = 0.0
    
    for f in fibers
        Cf_i = isotropic_stiffness(f[:E], f[:nu])
        S_i = eshelby_tensor_prolate(nu_m, f[:AR])
        
        # A_dil,i = [I + S_i * inv(Cm) * (Cf,i - Cm)]^-1
        Adil_i = inv(I6 + S_i * (inv(Cm) * (Cf_i - Cm)))
        
        push!(A_dil_list, Adil_i)
        v_total += f[:vf]
    end
    
    # 2. Calculate the common denominator for Mori-Tanaka
    # Denom = (1 - V_total)I + sum(V_j * Adil_j)
    sum_vA = zeros(6, 6)
    for (i, f) in enumerate(fibers)
        sum_vA += f[:vf] * A_dil_list[i]
    end
    MT_denom = inv((1.0 - v_total) * I6 + sum_vA)
    
    # 3. Sum up the contributions to the Aligned Stiffness
    C_aligned = Cm
    for (i, f) in enumerate(fibers)
        Cf_i = isotropic_stiffness(f[:E], f[:nu])
        A_MT_i = A_dil_list[i] * MT_denom
        C_aligned += f[:vf] * (Cf_i - Cm) * A_MT_i
    end
    
    # 4. Apply Orientation Averaging (reuse your existing Advani-Tucker logic)
    # ... call your existing averaging function using C_aligned ...
    
    return C_aligned
end


# # Define species with their respective volume fractions
# fibers = [
#     Dict(:name=>"Carbon", :vf=>0.21, :E=>230000.0, :nu=>0.2, :AR=>25.0),
#     Dict(:name=>"Glass",  :vf=>0.09, :E=>72000.0,  :nu=>0.22, :AR=>15.0)
#     ]

# # Run the hybrid homogenization
# C_hybrid_aligned = compute_hybrid_sfrp(2500.0, 0.35, fibers, 0.7, 0.2)

