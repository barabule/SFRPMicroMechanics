

##### FIBERPHASE


"""
    FiberPhase - collects fiber properties:
        - elastic_properties - (isotropic or Orthotropic) a subtype of AbstractElasticProperties
        - volume fraction of this phase
        - aspect ratio (ratiio of longest to shortest dim)
        - shape : InclusionGeometry (one of :SphericalInclusion(), SpheroidalInclusion(), NeedleInclusion(),
                    DiscInclusion(), ThinDiscInclusion()

"""
Base.@kwdef struct FiberPhase{EP<:AbstractElasticProperties, T<:Real, IG<:InclusionGeometry}   
    elastic_properties::EP
    volume_fraction::T
    aspect_ratio::T
    shape::IG = SpheroidalInclusion()

    function FiberPhase(ep, vf, ar, shap)
        VF, AR = promote(vf, ar)
        T = typeof(VF)
        EP = typeof(ep)
        IG = typeof(shap)
        @assert IG<:InclusionGeometry
        @assert EP<:AbstractElasticProperties
        @assert 0<= vf <= 1 "Volume fraction must be between 0 and 1!"
        @assert ar > 0 "Aspect ratio must be positive!"
        return new{EP, T, IG}(ep, VF, AR, shap)
    end
end


##### Constituent



abstract type Constituent end

@kwdef struct FiberConstituent{EP<:AbstractElasticProperties, 
                               T<:Real, 
                               OT<:AbstractOrientationTensor,
                               IG<:InclusionGeometry} <:Constituent
    elastic_properties<:EP
    density::T
    aspect_ratio::T
    shape::IG = SpheroidalInclusion()
    thermal_expansion<:ThermalExpansion{T}
end

#if only weight fraction is given, we cannot compute volume fraction without knowing the densities of fibers and matrix
#handle the orientation tensors and fraction independently from Constituents


##matrix

@kwdef struct MatrixConstituent{T<:Real} <: Constituent
    elastic_properties::IsotropicProperties{T}
    density::{T}
    thermal_expansion::ThermalExpansion{T}
end