# Isotropic Stiffness Matrix (Voigt 6x6)
function isotropic_stiffness(E, nu)
    λ = (E * nu) / ((1 + nu) * (1 - 2nu))
    μ = E / (2 * (1 + nu))

    return @SMatrix [λ+2μ   λ     λ     0   0   0;
                      λ     λ+2μ  λ     0   0   0;
                      λ     λ     λ+2μ  0   0   0;
                      0     0     0     μ   0   0;
                      0     0     0     0   μ   0;
                      0     0     0     0   0   μ]
end


"""
calc_vol_fraction(w_f, rho_f, rho_m)
Converts fiber weight fraction (e.g., 0.30 for 30% GF) to volume fraction.
"""
function calc_vol_fraction(w_f, rho_f, rho_m)
    v_f = (w_f / rho_f) / (w_f / rho_f + (1 - w_f) / rho_m)
    return v_f
end

function orthotropic_stiffness(; E1 = 1.0,
                             E2 = 1.0,
                             E3 = 1.0,
                             G12 = 1.0,
                            G23 = 1.0,
                            G31 = 1.0,
                            nu12 = nothing,
                            nu23 = nothing,
                            nu13 = nothing,
                            nu21 = nothing,
                            nu32 = nothing,
                            nu31 = nothing)

    @assert !(isnothing(nu12) && isnothing(nu21)) "ν12 or ν21 must be something"
    @assert !(isnothing(nu23) && isnothing(nu32)) "ν23 or ν32 must be something"
    @assert !(isnothing(nu13) && isnothing(nu31)) "ν31 or ν13 must be something"

    if isnothing(nu12)
        nu12 = nu21 * E1 / E2
    else
        nu21 = nu12 * E2 / E1
    end
    if isnothing(nu23)
        nu23 = nu32 * E2 / E3
    else
        nu32 = nu23 * E3 / E2
    end
    if isnothing(nu13)
        nu13 = nu31 * E1 / E3
    else
        nu31 = nu13 * E3 / E1
    end

    C = @SMatrix [    1/E1     -nu21/E2      -nu31/E3   0     0     0;
                    -nu12/E1     1/E2        -nu32/E3   0     0     0;
                    -nu13/E1   -nu23/E2        1/E3     0     0     0;
                       0          0             0     1/G23   0     0;
                       0          0             0       0    1/G31  0;
                       0          0             0       0     0    1/G12]

    # return(;stiffness = inv(C), 
    #         compliance = C)
    return C
    
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
    A_dil = inv(I6 + S * (inv(Cm) * (Cf - Cm)))
    A_MT = A_dil * inv((1 - vf) * I6 + vf * A_dil)
    return Cm + vf * (Cf - Cm) * A_MT
end

# function mori_tanaka(E_m, nu_m, E_f, nu_f, v_f, aspect_ratio)
#     Cm = isotropic_stiffness(E_m, nu_m)
#     Cf = isotropic_stiffness(E_f, nu_f)
#     I = Blueprints.eye(6) # Identity 6x6
    
#     # 1. Get Eshelby Tensor
#     S = eshelby_tensor_prolate(nu_m, aspect_ratio)
    
#     # 2. Dilute Concentration Tensor A_dil
#     # A_dil = [I + S * inv(Cm) * (Cf - Cm)]^-1
#     A_dil = inv(I + S * (inv(Cm) * (Cf - Cm)))
    
#     # 3. Mori-Tanaka Concentration Tensor A_MT
#     A_MT = A_dil * inv((1 - v_f) * I + v_f * A_dil)
    
#     # 4. Effective Stiffness
#     C_eff = Cm + v_f * (Cf - Cm) * A_MT
#     return C_eff
# end


"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
function hybrid_closure(orientation_tensor)
    a = Symmetric(SMatrix{3,3}(orientation_tensor))

    # 1. Calculate the interaction scalar 'f'
    # f = 1 - 27*det(a)
    # For a 3x3 orientation tensor, det(a) ranges from 0 (aligned) to 1/27 (random)
    f = 1.0 - 27.0 * det(a)
    
    # Ensure f stays within physical bounds due to numerical precision
    f = clamp(f, 0.0, 1.0)
    
    A_lin    = zeros(3, 3, 3, 3)
    A_quad   = zeros(3, 3, 3, 3)
    A_hybrid = zeros(3, 3, 3, 3)
    
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
    
    return A_hybrid
end


# Advani-Tucker Orientation Averaging
function orientation_average(C_aligned, a11, a22)
    a33 = 1.0 - a11 - a22
    a_mat = SMatrix{3,3}(diagm([a11, a22, a33]))
    A4 = hybrid_closure(a_mat)
    
    # Transversely Isotropic Invariants
    B1 = C_aligned[1,1] + C_aligned[2,2] - 2*C_aligned[1,2] - 4*C_aligned[6,6]
    B2 = C_aligned[1,2] - C_aligned[2,3]
    B3 = C_aligned[6,6] + 0.5*(C_aligned[2,3] - C_aligned[2,2])
    B4 = C_aligned[2,3]
    B5 = 0.5*(C_aligned[2,2] - C_aligned[2,3])
    
    C_avg = @MMatrix zeros(6,6)
    v = @SVector [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]

    δ(m,n) = (m==n) ? 1.0 : 0.0
    
    for r in 1:6, c in 1:6
        i,j = v[r]; k,l = v[c]
        C_avg[r,c] = B1*A4[i,j,k,l] + B2*(a_mat[i,j]*δ(k,l) + a_mat[k,l]*δ(i,j)) + 
                     B3*(a_mat[i,k]*δ(j,l) + a_mat[i,l]*δ(j,k) + a_mat[j,l]*δ(i,k) + a_mat[j,k]*δ(i,l)) + 
                     B4*(δ(i,j)*δ(k,l)) + B5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    end
    return C_avg
end


# Extract 9 Orthotropic Constants from Stiffness
function extract_orthotropic_constants(C_66)
    S = inv(C_66)
    E1, E2, E3 = 1/S[1,1], 1/S[2,2], 1/S[3,3]
    G23, G13, G12 = 1/S[4,4], 1/S[5,5], 1/S[6,6]
    return Dict("E1"=>E1, "E2"=>E2, "E3"=>E3, "G12"=>G12, "G23"=>G23, "G13"=>G13,
                "nu12"=>-S[2,1]*E1, "nu23"=>-S[3,2]*E2, "nu13"=>-S[3,1]*E1)
end


"""
    apparent_modulus(theta_deg, props)
    Calculates the Young's Modulus at an angle theta (degrees) 
    relative to the 1-axis.
    """
    function apparent_modulus(theta_deg, p)
        θ = deg2rad(theta_deg)
        s = sin(θ)
        c = cos(θ)
        
        # Off-axis compliance calculation
        inv_E_theta = (c^4 / p["E1"]) + 
            ( (1/p["G12"]) - (2*p["nu12"]/p["E1"]) ) * (s^2 * c^2) + 
                (s^4 / p["E2"])
            
            return 1.0 / inv_E_theta
    end
    
    
    """
    apparent_modulus_3d(phi_deg, theta_deg, C_avg_66)
    
    Calculates the Young's Modulus in 3D space.
    - phi_deg: Inclination from the 3-axis (0 to 180)
    - theta_deg: Azimuth from the 1-axis in the 1-2 plane (0 to 360)
    - C_avg_66: The 6x6 Voigt stiffness matrix from the module
    """
    function apparent_modulus_3d(phi_deg, theta_deg, C_avg_66)
        # 1. Convert C to 4th-order compliance S_ijkl
        S_66 = inv(C_avg_66)
        
        # Map Voigt 6x6 to 3x3x3x3 tensor S
        S = zeros(3, 3, 3, 3)
        v = [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]
        for r in 1:6, c in 1:6
            i, j = v[r]
            k, l = v[c]
            # Voigt shear components in S involve factors of 1, 1/2, or 1/4
            val = S_66[r, c]
            if r > 3; val /= 2.0; end
            if c > 3; val /= 2.0; end
            
            S[i,j,k,l] = S[j,i,k,l] = S[i,j,l,k] = S[j,i,l,k] = val
        end
        
        # 2. Define the direction unit vector n
        ϕ = deg2rad(phi_deg)
        θ = deg2rad(theta_deg)
        
        nx = sin(ϕ) * cos(θ)
        ny = sin(ϕ) * sin(θ)
        nz = cos(ϕ)
        n = [nx, ny, nz]
        
        # 3. Calculate 1/E = S_ijkl * ni * nj * nk * nl
        inv_E = 0.0
        for i=1:3, j=1:3, k=1:3, l=1:3
            inv_E += S[i,j,k,l] * n[i] * n[j] * n[k] * n[l]
        end
        
        return 1.0 / inv_E
    end


# # 3D Apparent Modulus Evaluation
# function apparent_modulus_3d(phi_deg, theta_deg, C_66)
#     S_66 = inv(C_66)
#     nx = sin(deg2rad(phi_deg)) * cos(deg2rad(theta_deg))
#     ny = sin(deg2rad(phi_deg)) * sin(deg2rad(theta_deg))
#     nz = cos(deg2rad(phi_deg))
#     n = [nx, ny, nz]
    
#     # 4th order compliance contraction
#     v = [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]
#     inv_E = 0.0
#     for r=1:6, c=1:6
#         i,j = v[r]; k,l = v[c]
#         val = S_66[r,c]
#         if r>3; val /= 2.0; end
#         if c>3; val /= 2.0; end
#         inv_E += val * n[i] * n[j] * n[k] * n[l]
#     end
#     return 1.0 / inv_E
# end

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
'fiber_list' is an array of dicts: [%, E, nu, AR]
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

"""
generate_stiffness_mesh(C_avg_66; resolution=50)

Generates a (nodes, faces) tuple for GLMakie.
- resolution: Number of points for theta and phi.
"""
function generate_stiffness_mesh(C_avg_66; resolution=60)
    # Define angles
    phi = range(0, stop=π, length=resolution)      # Inclination
    theta = range(0, stop=2π, length=resolution)    # Azimuth
    
    nodes = Point3f[]
    
    # 1. Generate Nodes
    # We use the previously defined 3D evaluation logic
    for p in phi
        for t in theta
            # Unit vector n
            nx = sin(p) * cos(t)
            ny = sin(p) * sin(t)
            nz = cos(p)
            
            # Evaluate E in this direction
            # Note: We call the apparent_modulus_3d logic here
            E_val = apparent_modulus_3d(rad2deg(p), rad2deg(t), C_avg_66)
            
            # Scale the unit vector by the modulus magnitude
            push!(nodes, Point3f(E_val * nx, E_val * ny, E_val * nz))
        end
    end
    
    # 2. Generate Faces (Triangulation of the sphere grid)
    faces = TriangleFace{Int}[]
    for i in 1:(resolution - 1)
        for j in 1:(resolution - 1)
            # Standard grid-to-mesh indexing
            p1 = (i - 1) * resolution + j
            p2 = i * resolution + j
            p3 = i * resolution + (j + 1)
            p4 = (i - 1) * resolution + (j + 1)
            
            push!(faces, TriangleFace(p1, p2, p3))
            push!(faces, TriangleFace(p1, p3, p4))
        end
    end
    
    return nodes, faces
end


"""
get_stiffness_slice(C_66, plane=:xy; resolution=100)

Returns (angles, moduli) for a specific plane.
    - :xy -> Plane 1-2 (Flow / Cross-flow)
    - :xz -> Plane 1-3 (Flow / Thickness)
    - :yz -> Plane 2-3 (Cross-flow / Thickness)
"""
function get_stiffness_slice(C_66, plane=:xy; resolution=180)
    alphas = range(0, 2π, length=resolution)
    moduli = Float64[]
    
    for α in alphas
        if plane == :xy
            # Inclination phi = 90 deg (XY plane), theta = alpha
            E = apparent_modulus_3d(90.0, rad2deg(α), C_66)
            elseif plane == :xz
            # theta = 0 deg, inclination phi = alpha
            E = apparent_modulus_3d(rad2deg(α), 0.0, C_66)
            elseif plane == :yz
            # theta = 90 deg, inclination phi = alpha
            E = apparent_modulus_3d(rad2deg(α), 90.0, C_66)
        else
            error("Plane must be :xy, :xz, or :yz")
        end
        push!(moduli, E)
    end
    return alphas, moduli
end