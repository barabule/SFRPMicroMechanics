###fiber

##### Constituent



abstract type Constituent end

@kwdef struct FiberConstituent{EP<:AbstractElasticProperties, 
                               T<:Real, 
                               OT<:AbstractOrientationTensor,
                               IG<:InclusionGeometry} <:Constituent
    elastic_properties::EP
    density::T
    aspect_ratio::T
    shape::IG = SpheroidalInclusion()
    thermal_expansion::ThermalExpansion{T}
end

#if only weight fraction is given, we cannot compute volume fraction without knowing the densities of fibers and matrix
#handle the orientation tensors and fraction independently from Constituents


##matrix

@kwdef struct MatrixConstituent{T<:Real} <: Constituent
    elastic_properties::IsotropicProperties{T}
    density::T
    thermal_expansion::ThermalExpansion{T}
end

## constructors 

function Constituent(;
                    elastic_properties = nothing::Union{Nothing, AbstractElasticProperties},
                    
                    density = nothing::Union{Nothing, Real},
                    aspect_ratio = nothing::Union{Nothing, Real},
                    shape = nothing::Union{Nothing, InclusionGeometry},
                    thermal_expansion = nothing::Union{Nothing, ThermalExpansion},
                    
                    )

    #if only isotropic elastic properties, density and thermal_expansion => MatrixConstituent
    if isa(elastic_properties, IsotropicProperties) && !isnothing(density) && 
                                                       !isnothing(thermal_expansion) && 
                                                       isnothing(orientation_tensor)
        # make thermal expansion isotropic
        cte = ThermalExpansion(thermal_expansion.alpha1)
        return MatrixConstituent(elastic_properties, density, thermal_expansion)
    end

    #only fiber
    @assert !isnothing(volume_fraction) && !isnothing(density)



end



"""
    effective_properties(matrix::MatrixConstituent, #matrix properties
                              fibers::AbstractVector{FiberConstituent}, #fibers 
                              fractions::AbstractVector{T}, #volume / weight fractions 
                              orientation_tensors::Union{AbstractOrientationTensor, AbstractVector{<:AbstractOrientationTensor}}; #...
                              by_weight = false, #fractions are per weight, default is by volume
                              mandel = true, #scale stiffness matrices with mandel or voigt factors
                              stiffness_method= :moritanaka, #method to compute stiffness
                              thermal_method = :moritanaka, #method to compute thermal expansion
                              closure_type = IBOF, #4th order orientation tensor approximation method
                              symmetrize = true,
                            ) where {T<:Real}

Computes effective properties of a composite:
    stiffness, thermal expansion and density
"""
function effective_properties(matrix::MatrixConstituent, #matrix properties
                              fibers::AbstractVector{FiberConstituent}, #fibers 
                              fractions::AbstractVector{T}, #volume / weight fractions 
                              orientation_tensors::Union{AbstractOrientationTensor, AbstractVector{<:AbstractOrientationTensor}}; #...
                              by_weight = false, #fractions are per weight, default is by volume
                              mandel = true, #scale stiffness matrices with mandel or voigt factors
                              stiffness_method= :moritanaka, #method to compute stiffness
                              thermal_method = :moritanaka, #method to compute thermal expansion
                              closure_type = IBOF, #4th order orientation tensor approximation method
                              symmetrize = true,
                            ) where {T<:Real}

    densities = vcat(matrix.density, [c.density for c in fibers])
    volume_fractions = by_weight ? to_volume_fractions(fractions, densities) : fractions

    if !isa(orientation_tensors, AbstractArray)
        avec = [orientation_tensors for _ in fibers]
    else
        @assert length(orientation_tensors) == length(fibers) "Length of fibers and a must the same!"
        avec = orientation_tensors
    end

    pm = matrix.elastic_properties
    fiber_phases = [FiberPhase(fibers[i].elastic_properties, 
                            volume_fractions[i], 
                            fibers[i].aspect_ratio,
                            fibers[i].shape) for i in eachindex(fibers)]

    #stiffness
    if stiffness_method == :moritanaka
        
        Ceff = mori_tanaka(pm, 
                    fiber_phases, 
                    orientation_tensors;
                    closure_type,
                    mandel,
                    symmetrize)

    elseif stiffness_method == :halpintsai
       Caligned = halpin_tsai(matrix, fibers,fractions; by_weight, mandel)
       #orientation_average
       a2_avg = isa(orientation_tensors, AbstractArray) ? average(orientation_tensors) : orientation_tensors
       Ceff = orientation_average(Caligned, a2_avg;closure_type, mandel)
    end
    peff = extract_orthotropic_constants(Ceff)
    
    #CTE
    ctes = vcat(matrix.thermal_expansion, [f.thermal_expansion for f in fibers])

    all_ctes = compute_all_thermal_expansions(pm, #matrix 
                                        fiber_phases,
                                        ctes,
                                        first(avec);
                                        mandel,
                                        average = true, #always orientation average
                                        )

    if thermal_method == :moritanaka
        cte_eff = all_ctes.moritanaka
    elseif thermal_method == :turner
        cte_eff = all.ctes.turner
    elseif thermal_method == :chow
        cte_eff = all_ctes.chow
    elseif thermal_method == :shapery
        cte_eff = all_ctes.shapery
    elseif thermal_method == :kerner
        cte_eff = all_ctes.kerner
    elseif thermal_method == :rom
        cte_eff = all_ctes.rom
    end

    # Density
    ρeff = effective_density(volume_fractions, densities)


    return (;density = ρeff,
            elastic_properties = peff,
            thermal_expansion = cte_eff,)


end