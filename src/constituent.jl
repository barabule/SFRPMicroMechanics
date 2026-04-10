###fiber

abstract type Constituent end

@kwdef struct FiberConstituent{EP<:AbstractElasticProperties, 
                               T<:Real, 
                               OT<:AbstractOrientationTensor} <:Constituent
    elastic_properties<:EP
    volume_fraction::T
    aspect_ratio::T
    shape<:InclusionGeometry
    density::T
    thermal_expansion<:ThermalExpansion{T}
    orientation_tensor::OT
end

#if only weight fraction is given, we cannot compute volume fraction without knowing the densities of fibers and matrix

##matrix

@kwdef struct MatrixConstituent{T<:Real} <: Constituent
    elastic_properties::IsotropicProperties{T}
    density::{T}
    thermal_expansion::ThermalExpansion{T}
end


## constructors 

function Constituent(;
                    elastic_properties = nothing::Union{Nothing, AbstractElasticProperties},
                    volume_fraction = nothing::Union{Nothing, Real},
                    weigth_fraction = nothing::Union{Nothing, Real},
                    density = nothing::Union{Nothing, Real},
                    aspect_ratio = nothing::Union{Nothing, Real},
                    shape = nothing::Union{Nothing, InclusionGeometry},
                    thermal_expansion = nothing::Union{Nothing, ThermalExpansion},
                    orientation_tensor = nothing::Union{Nothing, AbstractOrientationTensor},
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