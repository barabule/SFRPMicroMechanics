
# abstract type AbstractThermalExpansion end

struct ThermalExpansion{T<:Real}
    alpha1::T
    alpha2::T
    alpha3::T
    function ThermalExpansion(a1::T1, a2::T2,a3::T3) where {T1<:Real, T2<:Real, T3<:Real}
        args = promote(a1, a2, a3)
        T = eltype(args)
        return new{T}(args...)
    end
end

function ThermalExpansion(a1::T) where T<:Real
    return ThermalExpansion(a1, a1, a1)
end

function ThermalExpansion(a1::T1, a2::T2) where{T1<:Real, T2<:Real}
    return ThermalExpansion(a1, a2, a2)
end


function to_voigt(cte::ThermalExpansion{T}) where {T<:Real}
    return SVector{6, T}(cte.alpha1, cte.alpha2, cte.alpha3, 0, 0, 0)
end

function ThermalExpansion(v6::AbstractVector)
    @assert length(v6) == 6
    ThermalExpansion(v6[1], v6[2], v6[3])
end


function ThermalExpansion(matrix_properties, fiber_properties)
    """
    matrix = (;stiffness = Cm, 
                thermal_expansion = αm)
    fiber_properties = (;stiffness = [Cfi],
                        thermal_expansion = [αfi],
                        eshelby_tensors = [Sfi],
                        volume_fractions = [vfi])
    """


    Cm = matrix_properties.stiffness
    αm = matrix_properties.thermal_expansion |> to_voigt

    vf = fiber_properties.volume_fractions
    Sf = fiber_properties.eshelby_tensors
    Cf = fiber_properties.stiffness
    αf = fiber_properties.thermal_expansion

    I6 = SMatrix{6,6}(LinearAlgebra.I)

    term = sum((vf[i] * (Sf[i] - I6) for i in eachindex(vf)))
    
    αeff = MVector{6}(αm)
    for i in eachindex(vf)
        ΔC = Cf[i] - Cm
        αfi = to_voigt(αf[i])
        Δα = αfi - αm
        αeff .+= inv(ΔC * (Sf[i] - term) + Cm) * Cf[i] * Δα

    end
    return ThermalExpansion(αeff)
end

"""
    compute_sfrp_cte(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22)
Returns [alpha1, alpha2, alpha3] in the principal material directions.
"""
function ThermalExpansion(Em, num, alpham, Ef, nuf, alphaf, vf, AR, a11, a22) 
    # 1. Setup Stiffness
    Cm = isotropic_stiffness(Em, num)
    Cf = isotropic_stiffness(Ef, nuf)
    
    I6 = SMatrix{6, 6}(LinearAlgebra.I)

    inclusion = SpheroidalInclusion(num, AR)
    S_eshelby = eshelby_tensor(inclusion) |> convert_3333_to_66
    
    # 2. MT Concentration Tensor
    Adil = inv(I6 + S_eshelby * (inv(Cm) * (Cf - Cm)))
    Amt = Adil * inv((1 - vf) * I6 + vf * Adil)
    
    # 3. Aligned Effective Stiffness
    C_aligned = Cm + vf * (Cf - Cm) * Amt
    
    # 4. Aligned CTE (Rosen and Hashin logic)
    # Transformation to Voigt vector [a1, a2, a3, 0, 0, 0]
    alpha_m_vec = ThermalExpansion(alpham) |> to_voigt
    alpha_f_vec = ThermalExpansion(alphaf) |> to_voigt
    
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
    
    return ThermalExpansion(alpha1, alpha2, alpha3)
end