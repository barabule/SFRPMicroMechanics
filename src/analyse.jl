"""
    apparent_modulus(theta_deg, p)

apparent_modulus(theta_deg, p)
Calculates the Young's Modulus at an angle theta (degrees) 
relative to the 1-axis.
"""
function apparent_modulus(p::OrthotropicElasticParameters, theta_deg::T) where {T<:Real}
    

    # 1/E_theta = cos(theta)^4 / E1 + sin(theta)^4/E2 + (1/G12 - 2nu12/E1)*sin^2 cos^2

    s, c = sincosd(theta_deg)
    
    #nu12/E2 = nu21/E1
    E1 = p.E1
    E2 = p.E2
    nu12 = p.nu21 * E1 / E2
    G12 = p.G12
    # Off-axis compliance calculation
    inv_E_theta = (c^4 / E1) + 
                  (1 / G12 - 2*nu12 / E1) * (s^2 * c^2) + 
                  (s^4 / E2)
        
        return 1.0 / inv_E_theta
end
    
    
"""
    apparent_modulus(p::OrthotropicElasticParameters, theta_deg::T2, phi_deg::T1) where {T1<:Real, T2<:Real}

Calculates the Young's Modulus in 3D space.
- phi_deg: Inclination from the 3-axis (0 to 180)
- theta_deg: Azimuth from the 1-axis in the 1-2 plane (0 to 360)
- p - orthotropic elastic constants
"""
function apparent_modulus(p::OrthotropicElasticParameters, theta_deg::T1, phi_deg::T2) where {T1<:Real, T2<:Real}
    
    E1, E2, E3 = p.E1, p.E2, p.E3
    nu21, nu32, nu31 = p.nu21, p.nu32, p.nu31

    G12, G23, G31 = p.G12, p.G23, p.G31


    sϕ, cϕ = sincosd(phi_deg)
    sθ, cθ = sincosd(theta_deg)
    l1 = sϕ * cθ
    l2 = sϕ * sθ
    l3 = cϕ
    invE = l1^4 / E1 + l2^4 / E2 + l3^4 / E3 + 
           (1/G23 - 2nu32/E3) * l2^2 * l3^2 +
           (1/G31 - 2nu31/E3) * l1^2 * l3^2 +
           (1/G12 - 2nu21/E2) * l1^2 * l2^2

    
    return inv(invE)
end


function apparent_modulus(C66::AbstractMatrix, args...; kwargs...)
    return apparent_modulus(extract_orthotropic_constants(C66), args...; kwargs...)
end





"""
generate_stiffness_mesh(C_avg_66; resolution=50)

Generates a (nodes, faces) tuple for GLMakie.
- resolution: Number of points for theta and phi.
"""
function generate_stiffness_mesh(p::OrthotropicElasticParameters; resolution=60)
    C = stiffness_matrix_voigt(p)
    # Define angles
    phi = LinRange(0, 180, resolution)      # Inclination
    theta = LinRange(0, 360, resolution)    # Azimuth
    
    nodes = Point3f[]
    
    # 1. Generate Nodes
    # We use the previously defined 3D evaluation logic
    for p in phi
        for t in theta
            # Unit vector n
            sp, cp = sincosd(p)
            st, ct = sincosd(t)
            nx = sp * ct
            ny = sp * st
            nz = cp
            
            # Evaluate E in this direction
            # Note: We call the apparent_modulus_3d logic here
            E_val = apparent_modulus(C, t, p)
            
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