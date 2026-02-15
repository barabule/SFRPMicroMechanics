"""
    apparent_modulus(theta_deg, p)

apparent_modulus(theta_deg, p)
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

