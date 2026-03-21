

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


struct TransverseIsotropicElasticParameters{T<:Real}<:AbstractElasticParameters
    E1::T
    E2::T
    G12::T
    G23::T
    nu21::T
end

function stiffness_matrix_voigt(p::TransverseIsotropicElasticParameters; mandel = false)
    return stiffness_matrix_voigt(OrthotropicElasticParameters(p); mandel)    
end


function bulk_modulus(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    return E / (3(1-2nu))
end

function shear_modulus(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    return E/ (2(1+nu))
end

function lame_constant(p::IsotropicElasticParameters)
    return p.E_modulus * p.nu / ((1 + p.nu) * (1 - 2p.nu))
end

function OrthotropicElasticParameters(p::TransverseIsotropicElasticParameters)
    
    E1 = p.E1
    E2 = E3 = p.E2
    
    G31 = p.G12
    G23 = p.G23
    G12 = p.G12
    
    nu21 = nu31 = p.nu21
    nu23 = E2 / (2G23) - 1

    return OrthotropicElasticParameters(;E1, E2, E3, G12, G23, G31, nu21, nu31, nu23)
end


function OrthotropicElasticParameters(p::IsotropicElasticParameters)
    E, nu = p.E_modulus, p.nu
    G = E/(2(1+nu))
    return OrthotropicElasticParameters(;E1 = E,
                                         E2 = E,
                                         E3 = E,
                                         G12 = G,
                                         G23 = G,
                                         G31 = G,
                                         nu21 = nu,
                                         nu31 = nu,
                                         nu23 = nu)
end


"""
    compute_isotropic_best_fit(C_avg)
Projects a 6x6 stiffness tensor (Voigt notation) onto the isotropic manifold.
Returns IsotropicElasticParameters
"""
function isotropic_best_fit(C::AbstractMatrix)
    @assert size(C) == (6, 6)
    # Extract components for readability
    # Note: Using 1-based indexing for Julia
    C11, C22, C33 = C[1,1], C[2,2], C[3,3]
    C12, C23, C13 = C[1,2], C[2,3], C[1,3]
    C44, C55, C66 = C[4,4], C[5,5], C[6,6]

    # Best-fit Bulk Modulus (K)
    K = ((C11 + C22 + C33) + 2*(C12 + C23 + C13)) / 9.0
    
    # Best-fit Shear Modulus (G)
    G = ((C11 + C22 + C33) - (C12 + C23 + C13) + 3*(C44 + C55 + C66)) / 15.0

    # Convert to Young's Modulus (E) and Poisson's Ratio (nu)
    E = (9 * K * G) / (3 * K + G)
    nu = (3 * K - 2 * G) / (2 * (3 * K + G))

    return IsotropicElasticParameters(E, nu)
end

#upper bound
function voigt_average(ps::AbstractVector{<:AbstractElasticParameters}, w::AbstractVector{<:Real};mandel = true)
    @assert length(ps) == length(w)
    W =sum(w)
    return sum(stiffness_matrix_voigt(ps[i];mandel) * w[i]/ W for i in eachindex(ps))
end

#lower bound
function reuss_average(ps::AbstractVector{<:AbstractElasticParameters}, w::AbstractVector{<:Real};mandel = true)
    @assert length(ps) == length(w)
    W =sum(w)
    #maybe this could be be handled by compliance matrices
    return inv(sum(inv(stiffness_matrix_voigt(ps[i];mandel)) * w[i]/ W for i in eachindex(ps)))
end

#average upper lower
function hill_average(args...; kwargs...)
    1/2 * voigt_average(args...; kwargs...) + 1/2 * reuss_average(args...; kwargs...)
end


function IsotropicElasticParameters(ps::AbstractVector{<:AbstractElasticParameters}, w::AbstractVector{<:Real};mandel = true)
    @assert length(ps) == length(w)
    C_avg = hill_average(ps, w; mandel)
    return isotropic_best_fit(C_avg)
end



function Base.show(io::IO, ::MIME"text/plain", p::OrthotropicElasticParameters)
    println(io, "Elastic Constants:")
    # for field in fieldnames(T)
    #     value = getfield(p, field)
    #     println(io, "  $field = $value")
    # end
    println(io, "E1 = $(round(p.E1, sigdigits = 4))")
    println(io, "E2 = $(round(p.E2, sigdigits = 4))")
    println(io, "E3 = $(round(p.E3, sigdigits = 4))")

    println(io, "G12 = $(round(p.G12, sigdigits = 4))")
    println(io, "G31 = $(round(p.G31, sigdigits = 4))")
    println(io, "G23 = $(round(p.G23, sigdigits = 4))")


    nu12 = p.nu21 * p.E1 / p.E2
    nu13 = p.nu31 * p.E1 / p.E3
    nu23 = p.nu32 * p.E2 / p.E3

    println(io, "ν21 = $(round(p.nu21, sigdigits = 4))")
    println(io, "ν31 = $(round(p.nu31, sigdigits = 4))")
    println(io, "ν23 = $(round(nu23, sigdigits = 4))")
    println(io, "ν12 = $(round(nu12, sigdigits = 4))")
    println(io, "ν13 = $(round(nu13, sigdigits = 4))")
    println(io, "ν32 = $(round(p.nu32, sigdigits = 4))")

end

function Base.show(io::IO, ::MIME"text/plain", p::TransverseIsotropicElasticParameters)
    println(io, "Elastic Constants:")
    
    println(io, "E1 = $(round(p.E1, sigdigits = 4))")
    println(io, "E2 = $(round(p.E2, sigdigits = 4))")
    println(io, "E3 = $(round(p.E2, sigdigits = 4))")

    println(io, "G12 = $(round(p.G12, sigdigits = 4))")
    println(io, "G31 = $(round(p.G12, sigdigits = 4))")
    println(io, "G23 = $(round(p.G23, sigdigits = 4))")

    
    nu12 = p.nu21 * p.E1 / p.E2
    nu13 = nu12
    nu23 = p.E2 / (2p.G23) - 1
    nu32 = nu23

    println(io, "ν21 = $(round(p.nu21, sigdigits = 4))")
    println(io, "ν31 = $(round(p.nu21, sigdigits = 4))")
    println(io, "ν23 = $(round(nu23, sigdigits = 4))")
    println(io, "ν12 = $(round(nu12, sigdigits = 4))")
    println(io, "ν13 = $(round(nu13, sigdigits = 4))")
    println(io, "ν32 = $(round(p.nu32, sigdigits = 4))")
end

function Base.show(io::IO, ::MIME"text/plain", p::IsotropicElasticParameters)
    println(io, "Elastic Constants:")
    for field in fieldnames(T)
        value = getfield(p, field)
        println(io, "  $field = $(round(value, sigdigits =4))")
    end
end


function Base.isapprox(p1::IsotropicElasticParameters, p2::IsotropicElasticParameters; kwargs...)
    return isapprox(p1.E_modulus, p2.E_modulus; kwargs...) &&
           isapprox(p1.nu, p2.nu; kwargs...)
end


function Base.isapprox(p1::OrthotropicElasticParameters, p2::OrthotropicElasticParameters; kwargs...)
    return isapprox(p1.E1, p2.E1; kwargs...) &&
           isapprox(p1.E2, p2.E2; kwargs...) &&
           isapprox(p1.E3, p2.E3; kwargs...) &&
           isapprox(p1.G12, p2.G12; kwargs...) &&
           isapprox(p1.G23, p2.G23; kwargs...) &&
           isapprox(p1.G31, p2.G31; kwargs...) &&
           isapprox(p1.nu21, p2.nu21; kwargs...) &&
           isapprox(p1.nu32, p2.nu32; kwargs...) &&
           isapprox(p1.nu31, p2.nu31; kwargs...) 

end

function Base.isapprox(p1::TransverseIsotropicElasticParameters, p2::TransverseIsotropicElasticParameters; kwargs...)
    return isapprox(p1.E1, p2.E1; kwargs...) &&
           isapprox(p1.E2, p2.E2; kwargs...) &&
           isapprox(p1.G12, p2.G12; kwargs...) &&
           isapprox(p1.G23, p2.G23; kwargs...) &&
           isapprox(p1.nu21, p2.nu21; kwargs...)
end

function Base.isapprox(p1::OrthotropicElasticParameters, p2::AbstractElasticParameters; kwargs...) 
    return isapprox(p1, OrthotropicElasticParameters(p2);kwargs...)
end

Base.isapprox(p1::AbstractElasticParameters, p2::OrthotropicElasticParameters; kwargs...) = isapprox(p2,p1;kwargs...)

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


function TransverseIsotropicElasticParameters(;E1 = nothing,
                                               E2 = nothing,
                                               E3 = nothing,
                                               G12 = nothing,
                                               G23 = nothing,
                                               G31 = nothing,
                                               nu12 = nothing,
                                               nu21 = nothing,
                                               nu23 = nothing,
                                               nu32 = nothing,
                                               nu13 = nothing,
                                               nu31 = nothing)

    #we need only 5 constants but these have to describe the complete stiffness
    #longitudinal and transversal modulus
    @assert !isnothing(E1) "E1 must be given!"
    @assert !isnothing(E2) || !isnothing(E3) "Either E2 or E3 must be given"
    if isnothing(E2)
        E2 = E3
    end

    @assert !isnothing(G12) || !isnothing(G31) "Either G12 or G31 must be given!"
    if isnothing(G12)
        G12 = G31
    end
    
    @assert !isnothing(G23) || !isnothing(nu23) || !isnothing(nu32) "Either G23 or nu23 or nu32 must be given!"

    #in plane shear = isotropic
    if !isnothing(nu23)
        G23 = E2 / (2(1 + nu23))
    elseif !isnothing(nu32)
        G23 = E2 / (2(1 + nu32))
    end

    #LT plane 
    @assert !isnothing(nu12) || 
            !isnothing(nu21) || 
            !isnothing(nu13) || 
            !isnothing(nu31) "One of nu12, nu21, nu13 or nu31 must be given!"
    if !isnothing(nu12)
        nu21 = nu31 = nu12 * E2 / E1
    elseif !isnothing(nu13)
        nu21 = nu13 * E2 / E1
    elseif !isnothing(nu31)
        nu21 = nu31
    end


    nu12 = nu13 = nu21 * E1 / E2
    nu31 = nu21

    nu23 = nu32 = E2/(2G23) - 1
    Δ = 1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13
    # --- Stability Guards ---
    if Δ < 0
        error("Material is physically unstable (Δ = $Δ). Positive definiteness violated. Check Poisson's ratios.")
    elseif Δ < 1e-6
        @warn "Material is near a stability limit (Δ = $Δ). Results may be numerically sensitive."
    end
    return TransverseIsotropicElasticParameters(E1, E2, G12, G23, nu21)
end



function stiffness_matrix_voigt(p::IsotropicElasticParameters; mandel = false)
    E, nu = p.E_modulus, p.nu
    return isotropic_stiffness(E, nu; mandel)
end

function stiffness_matrix_voigt(p::OrthotropicElasticParameters; mandel = false)
    
    return orthotropic_stiffness(p.E1, p.E2, p.E3, p.G12, p.G23, p.G31, p.nu21, p.nu31, p.nu32; mandel) 
end


# Isotropic Stiffness Matrix (Voigt 6x6)
function isotropic_stiffness(E, nu; mandel = false)
    λ = (E * nu) / ((1 + nu) * (1 - 2nu))
    μ = E / (2 * (1 + nu))
    f = mandel ? 2 : 1
    # T = eltype((λ, f * μ))
    return SMatrix{6,6}(
                                [λ+2μ   λ     λ     0   0     0;
                                  λ     λ+2μ  λ     0   0     0;
                                  λ     λ     λ+2μ  0   0     0;
                                  0     0     0     f*μ 0     0;
                                  0     0     0     0   f*μ   0;
                                  0     0     0     0   0     f*μ])
end


# function stiffness_tensor(p::IsotropicElasticParameters{T}; mandel = false) where T
    
#     return mandel ? frommandel(stiffness_matrix_voigt(p;mandel)) : fromvoigt(stiffness_matrix_voigt(p;mandel))
# end


function is_structurally_isotropic(C66::AbstractMatrix)
    if size(C66) != (6,6) 
        return false
    end

    A = C66[1,1]
    B = C66[1,2]
    C = C66[4,4]
    M = @SMatrix  [A B B 0 0 0;
                   B A B 0 0 0;
                   B B A 0 0 0;
                   0 0 0 C 0 0;
                   0 0 0 0 C 0;
                   0 0 0 0 0 C]

    return all(M .≈ C66)
end


function is_isotropic(M::AbstractMatrix; tol=1e-9, mandel = false)
    size(M) == (6, 6) || throw(ArgumentError("Matrix must be 6x6"))
    
    # 1. Check for symmetry (Isotropic tensors are major-symmetric)
    if !isapprox(M, M', atol=tol)
        return false
    end

    # 2. Extract key components
    # Normal block (top-left 3x3)
    diag_normal = [M[1,1], M[2,2], M[3,3]]
    off_diag_normal = [M[1,2], M[1,3], M[2,3]]
    
    # Shear block (bottom-right 3x3)
    diag_shear = [M[4,4], M[5,5], M[6,6]]
    
    # 3. Check for uniformity within blocks
    # All diagonals of the normal block must be equal
    all_equal_diag = all(x -> isapprox(x, diag_normal[1], atol=tol), diag_normal)
    # All off-diagonals of the normal block must be equal
    all_equal_off = all(x -> isapprox(x, off_diag_normal[1], atol=tol), off_diag_normal)
    # All shear terms must be equal
    all_equal_shear = all(x -> isapprox(x, diag_shear[1], atol=tol), diag_shear)
    
    # 4. Check for zeros in coupling terms (off-diagonal blocks)
    # The 3x3 top-right and bottom-left blocks must be zero
    off_block_is_zero = isapprox(norm(M[1:3, 4:6]), 0, atol=tol) && 
                        isapprox(norm(M[4:6, 1:3]), 0, atol=tol)

    # 5. Check the Isotropy Relation
    # For Mandel: M44 = M11 - M12
    # For Voigt:  M44 = (M11 - M12) / 2
    rel_check = false
    if mandel
        rel_check = isapprox(diag_shear[1], (diag_normal[1] - off_diag_normal[1]), atol=tol)
    else
        rel_check = isapprox(diag_shear[1], (diag_normal[1] - off_diag_normal[1]) / 2, atol=tol)
    end

    return all_equal_diag && all_equal_off && all_equal_shear && off_block_is_zero && rel_check
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

function orthotropic_stiffness(E1, E2, E3, G12, G23, G31, nu21, nu31, nu32; mandel = false)
    nu12 = nu21 * E1 / E2
    nu13 = nu31 * E1 / E3
    nu23 = nu32 * E2 / E3
    f = mandel ? 2 : 1

    Δ = 1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13

    # --- Stability Guards ---
    if Δ < 0
        error("Material is physically unstable (Δ = $Δ). Positive definiteness violated. Check Poisson's ratios.")
    elseif Δ < 1e-6
        @warn "Material is near a stability limit (Δ = $Δ). Results may be numerically sensitive."
    end

    # 3. Stiffness Components
    # Formula check: C11 = E1(1 - nu23*nu32) / Δ
    C11 = E1 * (1 - nu23*nu32) / Δ
    C22 = E2 * (1 - nu31*nu13) / Δ
    C33 = E3 * (1 - nu12*nu21) / Δ

    C12 = E1 * (nu21 + nu31*nu23) / Δ
    C13 = E1 * (nu31 + nu21*nu32) / Δ
    C23 = E2 * (nu32 + nu12*nu31) / Δ

    # 4. Shear Components
    C44 = f * G23
    C55 = f * G31
    C66 = f * G12

    return @SMatrix [
        C11  C12  C13  0    0    0;
        C12  C22  C23  0    0    0;
        C13  C23  C33  0    0    0;
        0    0    0    C44  0    0;
        0    0    0    0    C55  0;
        0    0    0    0    0    C66
    ]  
end

function stiffness_tensor(p::OrthotropicElasticParameters;mandel=false)

    convert_66_to_3333(stiffness_matrix_voigt(p; mandel); mandel)
end


# Extract 9 Orthotropic Constants from Stiffness
function extract_orthotropic_constants(C_66::AbstractMatrix; mandel = true)
    @assert size(C_66) == (6,6) "Voigt matrix must be size 6x6!"
    S = inv(C_66) #compliance
    b = mandel ? 1/2 : 1
    E1, E2, E3 = 1/S[1,1], 1/S[2,2], 1/S[3,3]

    G23, G31, G12 = b/S[4,4], b/S[5,5], b/S[6,6]
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

function get_all_coefficients(a::OrientationTensor)
    a1 = a.a11
    a2 = a.a22
    a3 = 1 - a1 -a2
    a4 = 0
    a5 = 0
    a6 = 0
    return (a1, a2, a3, a4, a5, a6)
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

function get_all_coefficients(a::FullOrientationTensor)
    a3 = 1 - a.a1 - a.a2
    return (a.a1, a.a2, a3, a.a4, a.a5, a.a6)
end


function FullOrientationTensor(;a1 = nothing, a2 = nothing, a4 = nothing, a5 = nothing, a6 = nothing)
    @assert !(isnothing(a1) || isnothing(a2) || isnothing(a4) || isnothing(a5) || isnothing(a6))
    return FullOrientationTensor(a1, a2, a4, a5, a6)
end


function FullOrientationTensor(a::OrientationTensor{T}) where T
    a4, a5, a6 = 0, 0, 0
    return FullOrientationTensor(a.a11, a.a22, a4, a5, a6)
end


function decompose_eigenvalue(a::AbstractOrientationTensor)
    if isa(a, OrientationTensor)
        return (;tensor = a, rotation = one(SymmetricTensor{2,3}))
    end

    amat = to_matrix(a)
    lambda, vecs = eigen(amat) 
    # @info "lambda $lambda"
    idx = sortperm(lambda, rev= true)
    a11, a22 , _ = lambda[idx]
    
    R = Tensor{2,3}(vecs[:, idx])
    return (;tensor = OrientationTensor(a11, a22),
            rotation = R)
    
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



function to_matrix(a::AbstractOrientationTensor)
    a1, a2, a3, a4, a5,a6  = get_all_coefficients(a)
    return SymmetricTensor{2, 3}(
                    [a1 a6 a5;
                     a6 a2 a4;
                     a5 a4 a3])
end

function OrientationTensor(afull::FullOrientationTensor)
    lambda = LinearAlgebra.eigvals(to_matrix(afull))
    return OrientationTensor(lambda[3], lambda[2])
end


