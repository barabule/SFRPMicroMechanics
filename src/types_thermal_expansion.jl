

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