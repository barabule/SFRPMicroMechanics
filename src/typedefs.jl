

abstract type AbstractElasticParameters end

struct IsotropicElasticParameters{T<:Real}<:AbstractElasticParameters
    E_modulus::T
    nu::T
    function IsotropicElasticParameters(E, nu)
        args = promote(E, nu)
        T = eltype(args)
        @assert E>0 "E modulus must be positive!"
        return new{T}(args...)
    end
end


struct OrthotropicElasticParameters{T<:Real}<:AbstractElasticParameters
    E1::T
    E2::T
    E3::T
    G12::T
    G23::T
    G31::T
    nu21::T
    nu31::T
    nu32::T
    
    function OrthotropicElasticParameters(E1, E2, E3, G12, G23, G31, nu21, nu31, nu32)
            args = promote(E1, E2, E3, G12, G23, G31, nu21, nu31, nu32)
            T = eltype(args)
            
            @assert E1 > 0
            @assert E2 > 0
            @assert E3 > 0
            return new{T}(args...)
    end
end


function Base.show(io::IO, ::MIME"text/plain", p::T) where T<:AbstractElasticParameters
    println(io, "Elastic Constants:")
    for field in fieldnames(T)
        value = getfield(p, field)
        println(io, "  $field = $value")
    end
    if isa(p, OrthotropicElasticParameters)
        nu12 = p.nu21 * p.E1 / p.E2
        nu13 = p.nu31 * p.E1 / p.E3
        nu23 = p.nu32 * p.E2 / p.E3
        println("  nu12 = $nu12")
        println("  nu13 = $nu13")
        println("  nu23 = $nu23")
    end
end

function OrthotropicElasticParameters(;E1 = nothing,
                                       E2 = nothing,
                                       E3 = nothing,
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

    ## νij/Ei = νji / Ej
    if isnothing(nu21)
        nu21 = nu12 * E2 / E1
    end
    
    if isnothing(nu32)
        nu32 = nu23 * E3 / E2
    end

    if isnothing(nu31)
        nu31 = nu13 * E3 / E1
    end

    return OrthotropicElasticParameters(E1, E2, E3, G12, G23, G31, nu21, nu31, nu32)
end

function stiffness_matrix_voigt(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    return isotropic_stiffness(E, nu)
end

function stiffness_matrix_voigt(p::OrthotropicElasticParameters)
    
    return orthotropic_stiffness(p.E1, p.E2, p.E3, p.G12, p.G23, p.G31, p.nu21, p.nu31, p.nu32) 
end


# Isotropic Stiffness Matrix (Voigt 6x6)
function isotropic_stiffness(E, nu)
    λ = (E * nu) / ((1 + nu) * (1 - 2nu))
    μ = E / (2 * (1 + nu))

    return @SMatrix [λ+2μ   λ     λ     0   0   0;
                      λ     λ+2μ  λ     0   0   0;
                      λ     λ     λ+2μ  0   0   0;
                      0     0     0     μ   0   0;
                      0     0     0     0   μ   0;
                      0     0     0     0   0   μ]
end



function is_isotropic(C66)
    μ = C66[4,4]
    λ = C66[1,2]
    C66iso = @SMatrix [λ+2μ   λ     λ     0   0   0;
                      λ     λ+2μ  λ     0   0   0;
                      λ     λ     λ+2μ  0   0   0;
                      0     0     0     μ   0   0;
                      0     0     0     0   μ   0;
                      0     0     0     0   0   μ]
    return all(C66 .≈  C66iso)
end

"""
calc_vol_fraction(w_f, rho_f, rho_m)
Converts fiber weight fraction (e.g., 0.30 for 30% GF) to volume fraction.
"""
function calc_vol_fraction(mass_fraction, rho_f, rho_m)
    w_f = mass_fraction
    v_f = (w_f / rho_f) / (w_f / rho_f + (1 - w_f) / rho_m)
    return v_f
end

function orthotropic_stiffness(E1, E2, E3, G12, G23, G31, nu21, nu31, nu32)
    nu12 = nu21 * E1 / E2
    nu13 = nu31 * E1 / E3
    nu23 = nu32 * E2 / E3

    C = @SMatrix [    1/E1     -nu21/E2      -nu31/E3   0     0     0;
                    -nu12/E1     1/E2        -nu32/E3   0     0     0;
                    -nu13/E1   -nu23/E2        1/E3     0     0     0;
                       0          0             0     1/G23   0     0;
                       0          0             0       0    1/G31  0;
                       0          0             0       0     0    1/G12]

    return SMatrix(inv(C))    
end


# Extract 9 Orthotropic Constants from Stiffness
function extract_orthotropic_constants(C_66::AbstractMatrix)
    @assert size(C_66) == (6,6) "Voigt matrix must be size 6x6!"
    S = inv(C_66) #compliance

    E1, E2, E3 = 1/S[1,1], 1/S[2,2], 1/S[3,3]

    G23, G31, G12 = 1/S[4,4], 1/S[5,5], 1/S[6,6]
    nu12, nu23, nu13 = -S[2, 1] * E1, -S[3, 2] * E2, -S[3, 1] * E1
    
    return OrthotropicElasticParameters(;E1, E2, E3, G23, G31, G12, nu12, nu23, nu13)
   
end

abstract type AbstractOrientationTensor end

struct OrientationTensor{T<:Real} <:AbstractOrientationTensor
    a11::T
    a22::T

    function OrientationTensor(a11, a22)
        args = promote(a11, a22)
        T = eltype(args)
        @assert 1/3 <= a11 <= 1 "a11 must be between 1/3 and 1!"
        delta = sqrt(eps(T))
        @assert 1/2*(1 - a11)-delta <= a22 <= min(a11, 1-a11)+delta "a22 must be smaller than a11, and larger than a33!"
        return new{T}(args...)
    end
end

function to_matrix(A::OrientationTensor)
    a11, a22 = A.a11, A.a22
    a33 = 1 - a11 - a22
    return Symmetric(SMatrix{3,3}([a11 0 0;
                                    0 a22 0;
                                    0  0  a33]))
end


struct FullOrientationTensor{T<:Real} <:AbstractOrientationTensor
    a1::T
    a2::T #a3 is not needed because a1 + a2+ a3 =1
    a4::T
    a5::T
    a6::T

    function FullOrientationTensor(a1, a2, a4, a5, a6)
        args= promote(a1, a2, a4, a5, a6)
        T = eltype(args)
        @assert 1-a1-a2>=0 "a1 +a2 must be <=1 !"
        new{T}(args...)
    end
end

function FullOrientationTensor(;a1 = nothing, a2 = nothing, a4 = nothing, a5 = nothing, a6 = nothing)
    @assert !(isnothing(a1) || isnothing(a2) || isnothing(a4) || isnothing(a5) || isnothing(a6))
    return FullOrientationTensor(a1, a2, a4, a5, a6)
end


function Base.show(io::IO, ::MIME"text/plain", p::T) where {T<:AbstractOrientationTensor}
    println(io, "Orientation Tensor")
    if isa(p, OrientationTensor)
        println(io, "A1 = $(p.a11)")
        println(io, "A2 = $(p.a22)")
        println(io, "A2 = $(1 - p.a11 - p.a22)")
    elseif isa(p, FullOrientationTensor)
        println(io, "A1 = $(p.a1)")
        println(io, "A2 = $(p.a2)")
        println(io, "A3 = $(1 - p.a1 - p.a2)")
        println(io, "A4 = $(p.a4)")
        println(io, "A5 = $(p.a5)")
        println(io, "A6 = $(p.a6)")
    end
end


function to_matrix(a::FullOrientationTensor)
    return @SMatrix [a.a1 a.a6     a.a5;
                     a.a6 a.a2     a.a4;
                     a.a5 a.a4 (1- a.a1 - a.a2)]
end

function OrientationTensor(afull::FullOrientationTensor)
    lambda = LinearAlgebra.eigvals(to_matrix(afull))
    return OrientationTensor(lambda[3], lambda[2])
end
