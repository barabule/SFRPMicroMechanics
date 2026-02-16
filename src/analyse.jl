"""
    apparent_modulus(theta_deg, p)

apparent_modulus(theta_deg, p)
Calculates the Young's Modulus at an angle theta (degrees) 
relative to the 1-axis.
"""
function apparent_modulus(theta_deg, p::OrthotropicElasticParameters)
    θ = theta_deg

    # 1/E_theta = cos(theta)^4 / E1 + sin(theta)^4/E2 + (1/G12 - 2nu12/E1)*sin^2 cos^2

    s = sind(θ)
    c = cosd(θ)
    #nu12/E2 = nu21/E1
    E1 = p.E1
    E2 = p.E2
    nu12 = p.nu21 * E1 / E2
    G12 = p.G12
    # Off-axis compliance calculation
    inv_E_theta = (c^4 / E1) + 
        ( (1/G12) - (2*nu12/E1) ) * (s^2 * c^2) + 
            (s^4 / E2)
        
        return 1.0 / inv_E_theta
end
    
    
"""
    apparent_modulus_3d(phi_deg, theta_deg, p::OrthotropicElasticParameters)

Calculates the Young's Modulus in 3D space.
- phi_deg: Inclination from the 3-axis (0 to 180)
- theta_deg: Azimuth from the 1-axis in the 1-2 plane (0 to 360)
- p - orthotropic elastic constants
"""
function apparent_modulus_3d(phi_deg, theta_deg, p::OrthotropicElasticParameters)
    C_avg_66 = stiffness_matrix_voigt(p)
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


"""
generate_stiffness_mesh(C_avg_66; resolution=50)

Generates a (nodes, faces) tuple for GLMakie.
- resolution: Number of points for theta and phi.
"""
function generate_stiffness_mesh(p::OrthotropicElasticParameters; resolution=60)
    C_avg_66 = stiffness_matrix_voigt(p)
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


"""
    compute_sfrp_cte(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22)
Returns [alpha1, alpha2, alpha3] in the principal material directions.
"""
function compute_sfrp_cte(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22)
    # 1. Setup Stiffness
    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    I6 = Matrix{Float64}(I, 6, 6)
    S_eshelby = eshelby_tensor_spheroid(num, AR)
    
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

# #### MULTI FIBER
# """
# compute_hybrid_sfrp(matrix_props, fiber_list, a11, a22)
# 'fiber_list' is an array of dicts: [vf, E, nu, AR]
# """
# function compute_hybrid_sfrp(E_m, nu_m, fibers::Vector{AbstractElasticParameters}, A::OrientationTensor)
#     Cm = isotropic_stiffness(E_m, nu_m)
#     I6 = @SMatrix Matrix{Float64}(I, 6, 6)
    
#     # 1. Calculate Dilute Tensors for each fiber species
#     A_dil_list = []
#     v_total = 0.0
    
#     for f in fibers
#         Cf_i = stiffness_matrix_voigt(f)
#         S_i = eshelby_tensor_spheroid(nu_m, f[:AR])
        
        
#         Adil_i = inv(I6 + S_i * (inv(Cm) * (Cf_i - Cm)))
        
#         push!(A_dil_list, Adil_i)
#         v_total += f[:vf]
#     end
    
#     # 2. Calculate the common denominator for Mori-Tanaka
#     # Denom = (1 - V_total)I + sum(V_j * Adil_j)
#     sum_vA = zeros(6, 6)
#     for (i, f) in enumerate(fibers)
#         sum_vA += f[:vf] * A_dil_list[i]
#     end
#     MT_denom = inv((1.0 - v_total) * I6 + sum_vA)
    
#     # 3. Sum up the contributions to the Aligned Stiffness
#     C_aligned = Cm
#     for (i, f) in enumerate(fibers)
#         Cf_i = isotropic_stiffness(f[:E], f[:nu])
#         A_MT_i = A_dil_list[i] * MT_denom
#         C_aligned += f[:vf] * (Cf_i - Cm) * A_MT_i
#     end

#     return C_aligned
# end