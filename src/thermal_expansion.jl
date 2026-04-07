
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


function ThermalExpansion(pm::IsotropicElasticParameters, 
                          pf::IsotropicElasticParameters, 
                          cte_m::ThermalExpansion, 
                          cte_f::ThermalExpansion, 
                          vf::Real, 
                          AR::Real, 
                          a::AbstractOrientationTensor,
                          shape::InclusionGeometry;
                          mandel = false)
                          
    # Em, num = pm.E, pm.nu
    # Ef, nuf = pf.E, pf.nu

    alpham = cte_m.alpha1
    alphaf = cte_f.alpha1


    # 1. Setup Stiffness
    Cm = stiffness_matrix_voigt(pm; mandel)
    Cf = stiffness_matrix_voigt(pf; mandel)
    
    I6 = SMatrix{6, 6}(LinearAlgebra.I)

    
    S_eshelby = eshelby_tensor(shape, pm.nu, AR) |> convert_3333_to_66
    
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

    return orientation_averaging_thermal_expansion(ThermalExpansion(alpha_aligned_vec), a)
end


function orientation_averaging_thermal_expansion(cte::ThermalExpansion, a::OrientationTensor)


    αl = cte.alpha1
    αt = 1/2 * cte.alpha2 + 1/2 * cte.alpha3
    #averaging
    A2 = to_matrix(a)
    α_avg = SymmetricTensor{2,3}( 
                        (i, j) -> αt * δ(i, j) + (αl - αt) * A2[i, j]
                                )

    return  ThermalExpansion(tovoigt(α_avg))
end



function ThermalExpansion(pm::IsotropicElasticParameters, 
                          fiber::FiberPhase, 
                          ctes::Vector{<:ThermalExpansion}, 
                          a::AbstractOrientationTensor,
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



function ThermalExpansion(pm::IsotropicElasticParameters,
                            fibers::AbstractVector{<:FiberPhase},
                            ctes::AbstractVector{<:ThermalExpansion},
                            a::OrientationTensor;
                            mandel = true)

    @assert length(fibers) + 1 == length(ctes)
    return effective_thermal_expansion_mt(pm, fibers, ctes, a; mandel)

end


function effective_thermal_expansion_mt(pm::IsotropicElasticParameters, 
                                     fibers::AbstractVector{<:FiberPhase}, 
                                     ctes::AbstractVector{<:ThermalExpansion},
                                     a::AbstractOrientationTensor;
                                     mandel = true,
                                     average = true,
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

    #aligned CTE
    αeff = αm
    for i in 1:N
        αeff += ϕf[i] * inv((Cf[i] - Cm) * (Sf[i] - ΣϕfSI) + Cm) * Cf[i] * (αf[i] - αm)
    end

    cte_eff = ThermalExpansion(αeff)
    
    if average
        return orientation_averaging_thermal_expansion(cte_eff, a)
    end
    return cte_eff
end


"""
    effective_thermal_expansion_chow(pm::IsotropicElasticParameters, 
                                     pf::IsotropicElasticParameters, 
                                     cte_m::ThermalExpansion,
                                     cte_f::ThermalExpansion,
                                     volume_fraction::Real,
                                     aspect_ratio::Real,
                                     shape::InclusionGeometry)
                                     

Effective thermal expansion according to "Further studies on MorieTanaka models for thermal expansion
coefﬁcients of composites" by Pin Lu, 2012, Polymer 54, 1691-1699

Only works for isotropic fibers
"""
function effective_thermal_expansion_chow(pm::IsotropicElasticParameters, 
                                     pf::IsotropicElasticParameters, 
                                     cte_m::ThermalExpansion,
                                     cte_f::ThermalExpansion,
                                     volume_fraction::Real,
                                     aspect_ratio::Real,
                                     shape::InclusionGeometry;
                                     )

    
    Gm = shear_modulus(pm)
    Km = bulk_modulus(pm)

    
    
    Gf = shear_modulus(pf)
    Kf = bulk_modulus(pf)

    S = eshelby_tensor(shape, pm.nu, aspect_ratio)
    
    ϕf = volume_fraction
    
    αₘ = cte_m.alpha1
    γm = 3αₘ
    
    αf = cte_f
    γf = αf.alpha1 + αf.alpha2 + αf.alpha3

    b1 = S[2,2,2,2] + S[3,3,2,2] + S[1,1,2,2]
    b3 = S[2,2,1,1] + S[3,3,1,1] + S[1,1,1,1]

    c1 = S[2,2,2,2] + S[2,2,3,3] -2S[1,1,2,2]
    c3 = S[1,1,1,1] - S[2,2,1,1]

    K1 = 1 + (Kf/Km - 1) * ((1 - ϕf) * b1 + ϕf)
    K3 = 1 + (Kf/Km - 1) * ((1 - ϕf) * b3 + ϕf)

    G1 = 1 + (Gf/Gm - 1) * ((1 - ϕf) * c1 + ϕf)
    G3 = 1 + (Gf/Gm - 1) * ((1 - ϕf) * c3 + ϕf)

    α11 =       αₘ + ϕf * G1 / (2 * K1 * G3 + G1 * K3) * Kf / Km * (γf - γm)
    α22 = α33 = αₘ + ϕf * G3 / (2 * K1 * G3 + G1 * K3) * Kf / Km * (γf - γm)

    return ThermalExpansion(α11, α22, α33)
    
end




"""
    effective_thermal_expansion_shapery(pm::IsotropicElasticParameters, 
                                           pf::IsotropicElasticParameters, 
                                           cte_m::ThermalExpansion, 
                                           cte_f::ThermalExpansion, 
                                           volume_fraction::Real)
                                      
TBW
"""
function effective_thermal_expansion_shapery(pm::IsotropicElasticParameters, 
                                           pf::IsotropicElasticParameters, 
                                           cte_m::ThermalExpansion, 
                                           cte_f::ThermalExpansion, 
                                           volume_fraction::Real)

    Em, νm = pm.E, pm.nu
    Ef, νf = pf.E, pf.nu

    αm = cte_m.alpha1
    αf = cte_f.alpha1

    vf = volume_fraction
    vm = 1 - vf

    νc = vf * νf + vm * νf

    αL = (Ef * αf * vf + Em * αm * vm) / (Ef * vf + Em * vm)

    αT = (1 + νf) * αf * vf + (1 + νm) * αm * vm - αL * νc

    return ThermalExpansion(αL, αT)
end

function effective_thermal_expansion_kerner(pm::IsotropicElasticParameters,
                                            pf::IsotropicElasticParameters,
                                            cte_m::ThermalExpansion,
                                            cte_f::ThermalExpansion,
                                            vf::Real)

    Gm, Km = shear_modulus(pm), bulk_modulus(pm)
    Gf, Kf = shear_modulus(pf), bulk_modulus(pf)
    
    αm = cte_m.alpha1
    αf = cte_f.alpha1

    vm = 1 - vf

    αc = αm * vm + αf * vf + vm * (1 - vm) * (αf - αm) * (Kf - Km) / (vm * Km + vf * Kf + 3Km * Kf / (4Gm))

    return ThermalExpansion(αc)
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



function effective_thermal_expansion_mt_aligned(pm::IsotropicElasticParameters,
                                                pf::IsotropicElasticParameters,
                                                cte_m::ThermalExpansion,
                                                cte_f::ThermalExpansion,
                                                volume_fraction::Real,
                                                aspect_ratio::Real)

    λm = lame_constant(pm)
    μm = shear_modulus(pm)

    λf = lame_constant(pf)
    μf = shear_modulus(pf)

    S = eshelby_tensor(SpheroidalInclusion(), pm.nu, aspect_ratio)

    D1 = 1 + 2(μf - μm) / (λf - λm)
    D2 = (λm + 2μm) / (λf - λm)
    D3 = λm / (λf - λm)
    D0 = (3λf + 2μf) / (λf - λm)

    ϕf = volume_fraction
    
    B1 = ϕf * D1 + D2 + (1 - ϕf) * (D1 * S[1,1,1,1] + 2S[2,2,1,1])
    B2 = ϕf + D3 + (1 - ϕf) * (D1 * S[1,1,2,2] + S[2,2,2,2] + S[2,2,3,3])
    B3 = ϕf + D3 + (1 - ϕf) *(S[1,1,1,1] + (1 + D1) * S[2,2,1,1])
    B4 = ϕf * D1 + D2 + (1 - ϕf) * (S[1,1,2,2] + D1 * S[2,2,2,2] + S[2,2,3,3])
    B5 = ϕf + D3 + (1 - ϕf) * (S[1,1,2,2] + S[2,2,2,2] + D1 * S[2,2,3,3])

    αm = cte_m.alpha1
    αf = cte_f.alpha1 
    
    α11 =       αm + ϕf * (2B2 - B4 - B5) * D0 / (2B2*B3 - B1 * (B4 + B5)) * (αf - αm)
    α22 = α33 = αm + ϕf *       (B3 - B1) * D0 / (2B2*B3 - B1 * (B4 + B5)) * (αf - αm) 

    return ThermalExpansion(α11, α22, α33)
end


function effective_thermal_expansion_ROM(
                                        cte_m::ThermalExpansion,
                                        cte_f::ThermalExpansion,
                                        volume_fraction::Real,)

    αm = cte_m.alpha1
    αf = cte_f.alpha1

    ϕf = volume_fraction
    
    return ThermalExpansion(αm * (1 - ϕf) + αf * ϕf)
end


function effective_thermal_expansion_turner(pm::IsotropicElasticParameters,
                                        pf::IsotropicElasticParameters,
                                        cte_m::ThermalExpansion,
                                        cte_f::ThermalExpansion,
                                        volume_fraction::Real,)

    αm = cte_m.alpha1
    αf = cte_f.alpha1

    Km = bulk_modulus(pm)
    Kf = bulk_modulus(pf)

    ϕf = volume_fraction
    
    return ThermalExpansion((αm * (1 - ϕf) * Km + αf * ϕf * Kf) / ((1-ϕf) * Km + ϕf * Kf))
end



"""
    compute_all_thermal_expansions(pm::IsotropicElasticParameters, #matrix 
                                        fibers::AbstractVector{<:FiberPhase},
                                        ctes::AbstractVector{<:ThermalExpansion},
                                        a::AbstractOrientationTensor;
                                        mandel= true,
                                        average= false,
                                        )

Computes effective thermal expansion with several methods:
    ROM - rule of mixtures
    Turner
    Kerner
    Chow
    Shapery
    Moritanaka - 2 methods: multi fiber or simplified (1 fiber phase)

Usefull to have some idea of the range of plausible CTEs.


Returns a NamedTuple (;turner = cte_turner,
            rom = cte_rom,
            kerner = cte_kerner,
            shapery = cte_shapery,
            chow = cte_chow,
            moritanaka = cte_mt,
            moritanaka_simplified =cte_mt_simplified,
            )
                                        
"""
function compute_all_thermal_expansions(pm::IsotropicElasticParameters, #matrix 
                                        fibers::AbstractVector{<:FiberPhase},
                                        ctes::AbstractVector{<:ThermalExpansion},
                                        a::AbstractOrientationTensor;
                                        mandel= true,
                                        average= false,
                                        )

    @assert length(fibers) + 1 == length(ctes) "Length of matrix + fiber phases must be equal to length of ctes!"
    #if more than 1 fiber phase, compute an averaged "fiber" for some of the models
    cte_m = first(ctes)
    fiber_ctes = view(ctes,2:length(ctes))  
    f_avg = averaged_fiber(fibers, fiber_ctes)
    pf = f_avg.elastic_properties
    cte_f = f_avg.cte
    vf_avg = f_avg.volume_fraction
    AR_avg = f_avg.aspect_ratio


    #turner - 1 value
        cte_turner = effective_thermal_expansion_turner(pm, pf, cte_m, cte_f, vf_avg)
    #ROM 1 - value
        cte_rom = effective_thermal_expansion_ROM(cte_m, cte_f, vf_avg)
    #MT multi phase
        cte_mt = effective_thermal_expansion_mt(pm, fibers, ctes, a; mandel, average)
    #MT simplified
        cte_mt_simplified = effective_thermal_expansion_mt_aligned(pm, pf, cte_m, cte_f, vf_avg, AR_avg)
        if average
            cte_mt_simplified = orientation_averaging_thermal_expansion(cte_mt_simplified,a)
        end
    #kerner
        cte_kerner = effective_thermal_expansion_kerner(pm, pf, cte_m, cte_f, vf_avg)
    #shapery
        cte_shapery = effective_thermal_expansion_shapery(pm, pf, cte_m, cte_f, vf_avg)
        if average
            cte_shapery = orientation_averaging_thermal_expansion(cte_shapery, a)
        end
    #chow
        cte_chow = effective_thermal_expansion_chow(pm, pf, cte_m, cte_f, vf_avg, AR_avg, SpheroidalInclusion())
        if average
            cte_chow =orientation_averaging_thermal_expansion(cte_chow, a)
        end

    return (;turner = cte_turner,
            rom = cte_rom,
            kerner = cte_kerner,
            shapery = cte_shapery,
            chow = cte_chow,
            moritanaka = cte_mt,
            moritanaka_simplified =cte_mt_simplified,
            )
end




"""
    averaged_fiber(fibers::AbstractVector{<:FiberPhase}, 
                        ctes::AbstractVector{<:ThermalExpansion};
                        mandel = true)

Calculates a "best fit" isotropic fiber from the given fiber properties (elastic and thermal expansion)
Inputs:
    fibers - an AbstractVector of FiberPhase(s) containing elastic properties, volume fractions and aspect ratios
    ctes - an AbstractVector of ThermalExpansion(s)

kwargs: - mandel::Bool  - if true use Mandel scaling for offdiagonal elements, else Voigt

Returns a NamedTuple containing averaged (isotropic) properties:
    (;elastic_properties,
        cte::ThermalExpansion,
        volume_fraction,
        aspect_ratio)
"""
function averaged_fiber(fibers::AbstractVector{<:FiberPhase}, 
                        ctes::AbstractVector{<:ThermalExpansion};
                        mandel = true)

    @assert length(fibers) == length(ctes) "Length of fibers must be the same as ctes!"
    vfs = [fiber.volume_fraction for fiber in fibers]
    vf_tot = sum(vfs)
    @assert vf_tot <=1 "Sum of fiber volfractions must be <= 1!"

    αavg = 0.0
    for (cte, fiber) in zip(ctes, fibers)
        αavg += 1/3 * (cte.alpha1 + cte.alpha2 + cte.alpha3) * fiber.volume_fraction / vf_tot
    end
    
    cte_avg = ThermalExpansion(αavg)

    ps = [fiber.elastic_properties for fiber in fibers]
    w = vfs ./ vf_tot
    pf_avg = IsotropicElasticParameters(ps, w; mandel)

    AR_avg = sum([fiber.aspect_ratio for fiber in fibers] .* vfs ./ vf_tot)

    return (;elastic_properties = pf_avg,
             cte = cte_avg,
             volume_fraction = vf_tot,
             aspect_ratio = AR_avg)
end

