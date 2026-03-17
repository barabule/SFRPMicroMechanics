
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

function Base.show(io::IO, ::MIME"text/plain", p::T) where T<:ThermalExpansion
    println(io, "Thermal expansion coefficients:")
    println(io, "α₁ = $(p.alpha1)/°C")
    println(io, "α₂ = $(p.alpha2)/°C")
    println(io, "α₃ = $(p.alpha3)/°C")
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


# function ThermalExpansion(matrix_properties, fiber_properties)
#     """
#     matrix = (;stiffness = Cm, 
#                 thermal_expansion = αm)
#     fiber_properties = (;stiffness = [Cfi],
#                         thermal_expansion = [αfi],
#                         eshelby_tensors = [Sfi],
#                         volume_fractions = [vfi])
#     """


#     Cm = matrix_properties.stiffness
#     αm = matrix_properties.thermal_expansion |> to_voigt

#     vf = fiber_properties.volume_fractions
#     Sf = fiber_properties.eshelby_tensors
#     Cf = fiber_properties.stiffness
#     αf = fiber_properties.thermal_expansion

#     I6 = SMatrix{6,6}(LinearAlgebra.I)

#     term = sum((vf[i] * (Sf[i] - I6) for i in eachindex(vf)))
    
#     αeff = MVector{6}(αm)
#     for i in eachindex(vf)
#         ΔC = Cf[i] - Cm
#         αfi = to_voigt(αf[i])
#         Δα = αfi - αm
#         αeff .+= inv(ΔC * (Sf[i] - term) + Cm) * Cf[i] * Δα

#     end
#     return ThermalExpansion(αeff)
# end


function ThermalExpansion(pm::IsotropicElasticParameters, 
                          pf::IsotropicElasticParameters, 
                          cte_m::ThermalExpansion, 
                          cte_f::ThermalExpansion, 
                          vf::Real, 
                          AR::Real, 
                          a::OrientationTensor,
                          shape::InclusionGeometry;
                          mandel = false)
                          
    Em, num = pm.E_modulus, pm.nu
    Ef, nuf = pf.E_modulus, pf.nu

    alpham = cte_m.alpha1
    alphaf = cte_f.alpha1

    a11, a22 = a.a11, a.a22

    # 1. Setup Stiffness
    Cm = stiffness_matrix_voigt(pm; mandel)
    Cf = stiffness_matrix_voigt(pf; mandel)
    
    I6 = SMatrix{6, 6}(LinearAlgebra.I)

    
    S_eshelby = eshelby_tensor(shape, num, AR) |> convert_3333_to_66
    
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


function ThermalExpansion(pm::IsotropicElasticParameters, 
                          fiber::FiberPhase, 
                          ctes::Vector{<:ThermalExpansion}, 
                          a::OrientationTensor,
                          )

    @assert length(ctes)>=2 "2 CTEs must be given!"
    @assert isa(fiber.elastic_properties, IsotropicElasticParameters) "Fiber must be isotropic!"

    return ThermalExpansion(pm, 
                            fiber.elastic_properties,
                            ctes[1],
                            ctes[2],
                            fiber.volume_fraction,
                            fiber.aspect_ratio,
                            a,
                            fiber.shape)

end


# function ThermalExpansion(pm::IsotropicElasticParameters, 
#                           pf::TransverseIsotropicElasticParameters,
#                           cte_m::ThermalExpansion, 
#                           cte_f::ThermalExpansion, 
#                           vf::Real, 
#                           AR::Real, 
#                           a::OrientationTensor,
#                           shape::InclusionGeometry;
#                           mandel = true)
    
#     # 1. Setup Phase Properties
#     Cm = stiffness_matrix_voigt(pm;mandel)
#     Cf = stiffness_matrix_voigt(pf;mandel) 
    
#     # 2. Micromechanics (Mori-Tanaka)
#     I6 = SMatrix{6, 6}(LinearAlgebra.I)
#     S  = eshelby_tensor(shape, pm.nu, AR) |> convert_3333_to_66
    
#     # Concentration tensor (Dilute -> Mori-Tanaka)
#     A_dil = inv(I6 + S * (inv(Cm) * (Cf - Cm)))
#     A_mt  = A_dil * inv((1 - vf) * I6 + vf * A_dil)
    
#     # Aligned Effective Stiffness
#     C_aligned = Cm + vf * (Cf - Cm) * A_mt

#     # 3. Levin's Relation for Aligned CTE Vector
#     # α_aligned = α_m + vf * [inv(C_aligned - Cm) * (Cf - Cm) * A_mt] * (α_f - α_m)
#     am_vec = to_voigt(cte_m)
#     af_vec = to_voigt(cte_f)
    
#     term = inv(C_aligned - Cm) * (Cf - Cm) * A_mt * (af_vec - am_vec)
#     alpha_aligned = am_vec + vf * term

#     # 4. Orientation Averaging (2nd Order Tensor)
#     # Since fiber is Transversely Isotropic: α_long = α[1], α_trans = α[2] = α[3]
#     α_l, α_t = alpha_aligned[1], alpha_aligned[2]
    
#     a11, a22 = a.a11, a.a22
#     a33 = 1.0 - a11 - a22
    
#     # Map back to the global coordinate system using orientation components
#     return ThermalExpansion(
#         (α_l - α_t) * a11 + α_t,
#         (α_l - α_t) * a22 + α_t,
#         (α_l - α_t) * a33 + α_t
#     )
# end

function ThermalExpansion(pm::IsotropicElasticParameters,
                            fibers::AbstractVector{<:FiberPhase},
                            ctes::AbstractVector{<:ThermalExpansion};
                            mandel = true)

    @assert length(fibers) + 1 == length(ctes)
    return effective_thermal_expansion_mt(pm, fibers, ctes; mandel)

end


function effective_thermal_expansion_mt(pm::IsotropicElasticParameters, 
                                     fibers::AbstractVector{<:FiberPhase}, 
                                     ctes::AbstractVector{<:ThermalExpansion};
                                     mandel = true,
                                     )

    N = length(fibers)
    @assert N+1 == length(ctes)

    αm = to_voigt(first(ctes))
    αf = [to_voigt(ctes[i]) for i in 2:lastindex(ctes)]

    νm = pm.nu

    Cm = stiffness_matrix_voigt(pm; mandel)
    T = eltype(Cm)
    ϕf = [f.volume_fraction for f in fibers]
    Cf = [stiffness_matrix_voigt(f.elastic_properties; mandel) for f in fibers]
    Sf = [convert_3333_to_66(eshelby_tensor(f.shape, νm, f.aspect_ratio);mandel) for f in fibers]


    ΣϕfSI = SMatrix{6,6}(zeros(T, 6,6))
    for i in 1:N
        ΣϕfSI += ϕf[i] * (Sf[i] - LinearAlgebra.I)
    end
    # Σϕf = sum(ϕf)


    αeff = αm
    for i in 1:N
        αeff += ϕf[i] * inv((Cf[i] - Cm) * (Sf[i] - ΣϕfSI) + Cm) * Cf[i] * (αf[i] - αm)
    end
    # @info αeff
    return ThermalExpansion(αeff)
end


function effective_thermal_expansion_chow(pm::IsotropicElasticParameters, 
                                     fiber::FiberPhase, 
                                     ctes::AbstractVector{<:ThermalExpansion})

    νm = pm.nu
    Em = pm.E_modulus
    Gm = Em / (2 * (1 + νm))
    Km = Em / (3 * (1- 2νm))
    pf = fiber.elastic_properties
    νf = pf.nu
    Ef = pf.E_modulus
    Gf = Ef/(2 * (1+νf))
    Kf = Ef / (3 * (1- 2νf))

    S = eshelby_tensor(fiber.shape, νf, fiber.aspect_ratio)
    ϕf = fiber.volume_fraction
    αₘ = ctes[1].alpha1
    γm = 3αₘ
    αf = ctes[2]
    γf = αf.alpha1 + αf.alpha2 + αf.alpha3


    b1 = S[2,2,2,2] + S[3,3,2,2] + S[1,1,2,2]
    b3 = S[2,2,1,1] + S[3,3,1,1] + S[1,1,1,1]

    c1 = S[2,2,2,2] + S[2,2,3,3] -2S[1,1,2,2]
    c3 = S[1,1,1,1] - S[2,2,1,1]

    K1 = 1 + (Kf/Km - 1) * ((1 - ϕf) * b1 + ϕf)
    K3 = 1 + (Kf/Km - 1) * ((1 - ϕf) * b3 + ϕf)

    G1 = 1 + (Gf/Gm - 1) * ((1 - ϕf) * c1 + ϕf)
    G3 = 1 + (Gf/Gm - 1) * ((1 - ϕf) * c3 + ϕf)

    α11 = αₘ + ϕf * G1 / (2 * K1 * G3 + G1 * K3 * Km) * Kf / Km * (γf - γm)
    α22 = α33 = αₘ + ϕf * G3 / (2 * K1 * G3 + G1 * K3 * Km) * Kf / Km * (γf - γm)

    return ThermalExpansion(α11, α22, α33)
end




function compute_thermal_expansion_shapery(pm::IsotropicElasticParameters, 
                                           pf::IsotropicElasticParameters, 
                                           cte_m::ThermalExpansion, 
                                           cte_f::ThermalExpansion, 
                                           vf::Real)

    Em, νm = pm.E_modulus, pm.nu
    Ef, νf = pf.E_modulus, pf.nu

    αm = cte_m.alpha1
    αf = cte_f.alpha1

    vm = 1 - vf

    νc = vf * νf + vm * νf

    αL = (Ef * αf * vf + Em * αm * vm) / (Ef * vf + Em * vm)

    αT = (1 + νf) * αf * vf + (1 + νm) * αm * vm - αL * νc

    return ThermalExpansion(αL, αT, αT)
end



function compute_thermal_expansion_kerner(pm::IsotropicElasticParameters, 
                                          pfs::Vector{IsotropicElasticParameters}, 
                                          cte_m::ThermalExpansion, 
                                          ctes_f::Vector{ThermalExpansion}, 
                                          vfs::Vector{Real})

        αm = cte_m.alpha1
        Vtol = sum(vfs)
        Vm = 1 - Vtol

        αp = sum(cts_f[i].alpha1 * vfs[i]/Vtol for i in eachindex(ctes_f))

        Km = bulk_modulus(pm)
        Gm = shear_modulus(pm)
        Kp = sum(bulk_modulus.(pfs) .* vfs ./ Vtol)
        # Gp = sum(shear_modulus(pfs[i]) * vfs[i] / Vtol for i in eachindex(pfs))

        αC = αm * Vm + sum(ctes_f[i].alpha1 * vfs[i] for i in eachindex(ctes_f)) + Vm * (1 - Vm) * 
            (αp - αm) * (Kp - Km) / (Vm * Km + sum(vfs[i] .* Kp) + (3 * Km * Kp) / (4Gm))

        return ThermalExpansion(αC)

end

function bulk_modulus(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    return E / (3(1-2nu))
end

function shear_modulus(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    return E/ (2(1+nu))
end


function hashin_strikman_bounds(pm::IsotropicElasticParameters, pf::IsotropicElasticParameters, vf)
    vm = 1 - vf

    Km = bulk_modulus(pm)
    Gm = shear_modulus(pm)
    Kf = bulk_modulus(pf)
    Gf = shear_modulus(pf)

    Kl = Km + vf/ (1/(Kf - Km) + vm / (Km + 4Gm/3))

    Ku = Kf + vm / (1 / (Km - Kf) + (1 - vm) / (Kf + 4Gf/3))

    return (;lower= Kl, upper = Ku)
end

