module MoriTanaka

using LinearAlgebra

"""
    stiffness_isotropic(E, nu)
Returns the 6x6 Voigt stiffness matrix for an isotropic material.
Order: [11, 22, 33, 23, 13, 12]
"""
function stiffness_isotropic(E, nu)
    λ = (E * nu) / ((1 + nu) * (1 - 2nu))
    μ = E / (2 * (1 + nu))
    
    C = zeros(6, 6)
    # Longitudinal/Transverse blocks
    C[1,1] = C[2,2] = C[3,3] = λ + 2μ
    C[1,2] = C[1,3] = C[2,1] = C[2,3] = C[3,1] = C[3,2] = λ
    # Shear components
    C[4,4] = C[5,5] = C[6,6] = μ
    return C
end

"""
    eshelby_spheroid(ar, nu)
Calculates the Eshelby tensor S for a prolate spheroid (aspect ratio ar > 1)
aligned along the 1-axis in an isotropic matrix with Poisson's ratio nu.
"""
function eshelby_spheroid(ar, nu)
    # Auxiliary variable g for prolate spheroids
    # ar = length / diameter
    g = (ar / (ar^2 - 1)^1.5) * (ar * sqrt(ar^2 - 1) - acosh(ar))
    
    # Components S_ijkl in Voigt notation
    S = zeros(6, 6)
    
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
    
    return S
end

"""
    mori_tanaka_homogenization(Cm, Ci, S, f)
Cm: Matrix Stiffness, Ci: Inclusion Stiffness, S: Eshelby, f: Volume fraction
"""
function mori_tanaka_homogenization(Cm, Ci, S, f)
    I = diagm(ones(6))
    # Dilute concentration tensor (T)
    T = inv(I + S * (inv(Cm) * (Ci - Cm)))
    
    # Mori-Tanaka concentration tensor (A)
    A = T * inv((1 - f) * I + f * T)
    
    # Effective Stiffness
    C_eff = Cm + f * (Ci - Cm) * A
    return C_eff
end

# Extract 9 Orthotropic Constants from Stiffness
function extract_orthotropic_constants(C_66::AbstractMatrix)
    @assert size(C_66) == (6,6) "Voigt matrix must be size 6x6!"
    S = inv(C_66) #compliance

    E1, E2, E3 = 1/S[1,1], 1/S[2,2], 1/S[3,3]

    G23, G31, G12 = 1/S[4,4], 1/S[5,5], 1/S[6,6]
    nu12, nu23, nu13 = -S[2, 1] * E1, -S[3, 2] * E2, -S[3, 1] * E1
    return (;E1, E2, E3, G23, G31, G12, nu12, nu23, nu13)
   
end


end #module


MT = MoriTanaka

function run_me(; E_m= 3.0,
                    nu_m = 0.35,
                    E_f = 72.0,
                    nu_f = 0.22,
                    vf = 0.15,
                    AR = 20.0,
                    )
    # # --- RUN EXAMPLE ---
    # # Polyamide 6 (Matrix) + Glass Fiber (Inclusion)
    # E_m, nu_m = 3.0, 0.35    # GPa
    # E_f, nu_f = 72.0, 0.22   # GPa
    # vol_fraction = 0.15        # 15% fiber
    # aspect_ratio = 25.0        # Short fibers

    Cm = MT.stiffness_isotropic(E_m, nu_m)
    Ci = MT.stiffness_isotropic(E_f, nu_f)
    S  = MT.eshelby_spheroid(AR, nu_m)

    C_star = MT.mori_tanaka_homogenization(Cm, Ci, S, vf)

    println("Homogenized Stiffness Matrix (GPa):")
    display(round.(C_star, digits=2))
    display(MT.extract_orthotropic_constants(C_star))
end
