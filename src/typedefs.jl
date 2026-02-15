

abstract type AbstractElasticParameters end

struct IsotropicElasticParameters{T<:Real}
    E_modulus::T
    nu::T
    function IsotropicElasticParameters(E, nu)
        args = promote(E, nu)
        T = eltype(args)
        @assert E>0 "E modulus must be positive!"
        return new{T}(args...)
    end
end


struct OrthotropicElasticParameters{T<:Real}
    E1::T
    E2::T
    E3::T
    G12::T
    G23::T
    G31::T
    nu12::T
    nu21::T
    nu23::T
    nu32::T
    nu13::T
    nu31::T

    function OrthotropicElasticParameters(E1, E2, E3, G12, G23, G31, nu12, nu21, nu23, nu32, nu13, nu31)
            args = promote(E1, E2, E3, G12, G23, G31, nu12, nu21, nu23, nu32, nu13, nu31)
            T = eltype(args)
            @assert nu12 / E1 ≈ nu21 / E2
            @assert nu23 / E2 ≈ nu32 / E3
            @assert nu13 / E1 ≈ nu31 / E3
            @assert E1 > 0
            @assert E2 > 0
            @assert E3 > 0
            return new{T}(args...)
    end
end


function OrthotropicElasticParameters(;E1=nothing,
                                       E2=nothing,
                                       E3=nothing,
                                       nu12 = nothing,
                                       nu21 = nothing,
                                       nu13 = nothing,
                                       nu31 = nothing,
                                       nu23 = nothing,
                                       nu32 = nothing,
                                       G12 = nothing,
                                       G23 = nothing,
                                       G31 = nothing,
                                       )

    if isnothing(E1) || isnothing(E2) || isnothing(E3) || isnothing(G12) || isnothing(G23) || isnothing(G31)
        error("Elastic moduli must not be nothing!")
    end

    @assert !(isnothing(nu12) && isnothing(nu21)) "Either nu12 or nu21 must be given!"
    @assert !(isnothing(nu23) && isnothing(nu32)) "Either nu23 or nu32 must be given!"
    @assert !(isnothing(nu31) && isnothing(nu13)) "Either nu31 or nu13 must be given!"

    if isnothing(nu12)
        nu12 = nu21 * E1 / E2
    else
        nu21 = nu12 * E2 / E1
    end

    if isnothing(nu23)
        nu23 = nu32 * E2 / E3
    else
        nu32 = nu23 * E3 / E2
    end

    if isnothing(nu13)
        nu13 = nu31 * E1 / E3
    else
        nu31 = nu13 * E3 / E1
    end

    return OrthotropicElasticParameters(E1, E2, E3, G12, G23, G31, nu12, nu21, nu23, nu32, nu13, nu31)
end

function stiffness_matrix_voigt(elastic_parameters::IsotropicElasticParameters)
    E, nu = elastic_parameters
    return isotropic_stiffness(E, nu)
end

function stiffness_matrix_voigt(elastic_parameters::OrthotropicElasticParameters)
    E1, E2, E3, G12, G23, G31, nu12, nu21, nu23, nu32, nu13, nu31 = elastic_parameters
    return orthotropic_stiffness(E1, E2, E3, G12, G23, G31, nu12, nu21, nu23, nu32, nu13, nu31) 
end


struct OrientationTensor{T<:Real}
    a11::T
    a22::T

    function OrientationTensor(a11, a22)
        T = promote_type(a11, a22)
        @assert 1/3 <= a11 <= 1 "a11 must be between 1/3 and 1!"
        @assert 1/2*(1 - a11) <= a22 <= min(a11, 1-a11) "a22 must be smaller than a11, and larger than a33!"
        return new{T}(a11, a22)
    end
end