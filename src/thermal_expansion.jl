
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