"""
HomoPy.jl — Julia translation of the HomoPy Python package by Nicolas Christ et al.
Original: https://github.com/Extraweich/homopy
License: MIT

Implements:
  1. Halpin-Tsai (with Shear-Lag / Cox modification) — scalar homogenization
  2. Mori-Tanaka (Benveniste formulation + optional Segura symmetrization) — tensorial homogenization

Dependencies: LinearAlgebra, StaticArrays
"""

module HomoPy

using LinearAlgebra
using StaticArrays

export isotropic_stiffness, transverse_isotropic_stiffness
export rotate_stiffness_4th, mandel_to_4th, fourth_to_mandel
export eshelby_tensor_needle, eshelby_tensor_spheroid, eshelby_tensor_sphere
export mori_tanaka, mori_tanaka_multi, mori_tanaka_orientation
export halpin_tsai, halpin_tsai_UD
export young_modulus_from_stiffness, compliance_from_stiffness

# ──────────────────────────────────────────────────────────────────────────────
# Constants / index helpers
# ──────────────────────────────────────────────────────────────────────────────

"""kron delta δ(i, j) = 1 if i==j else 0"""
@inline δ(i, j) = i == j 


"""Voigt notation mapping (i,j) -> Voigt index (1-based)."""
const VOIGT = SMatrix{3,3,Int}([1 6 5; 
                                6 2 4; 
                                5 4 3])

"""Mandel scaling factors for off-diagonal entries."""
mandel_scale(i, j) = δ(i, j) ? 1.0 : sqrt(2.0)

# ──────────────────────────────────────────────────────────────────────────────
# Stiffness tensor constructors
# ──────────────────────────────────────────────────────────────────────────────

"""
    isotropic_stiffness(E, ν) -> SMatrix{6,6}

Build 6×6 Voigt stiffness matrix for an isotropic material from
Young's modulus `E` and Poisson ratio `ν`.
"""
function isotropic_stiffness(E::Real, ν::Real)
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    C = @MMatrix zeros(6, 6)
    for i in 1:3, j in 1:3
        C[i, j] = λ
    end
    for i in 1:3
        C[i, i] += 2μ
    end
    for i in 4:6
        C[i, i] = μ
    end
    return SMatrix{6,6}(C)
end

"""
    transverse_isotropic_stiffness(E_L, E_T, ν_L, ν_TT, G_L) -> SMatrix{6,6}

Build 6×6 Voigt stiffness matrix for a transversely isotropic material.
The longitudinal axis is x₁ (fiber direction).

Arguments:
  E_L   — longitudinal Young's modulus
  E_T   — transverse Young's modulus
  ν_L   — longitudinal Poisson ratio (ν₁₂ = ν₁₃)
  ν_TT  — transverse Poisson ratio (ν₂₃)
  G_L   — longitudinal shear modulus (G₁₂ = G₁₃)
"""
function transverse_isotropic_stiffness(E_L::Real, E_T::Real, ν_L::Real, ν_TT::Real, G_L::Real)
    G_TT = E_T / (2(1 + ν_TT))
    Δ = (1 - ν_TT^2 * E_L / E_T - 2 * ν_L^2 * E_T / E_L - 2 * ν_L^2 * ν_TT * E_T / E_L) / (E_L * E_T^2)
    # compliance
    S = @MMatrix zeros(6, 6)
    S[1, 1] = 1 / E_L
    S[2, 2] = 1 / E_T
    S[3, 3] = 1 / E_T
    S[1, 2] = S[2, 1] = -ν_L / E_L
    S[1, 3] = S[3, 1] = -ν_L / E_L
    S[2, 3] = S[3, 2] = -ν_TT / E_T
    S[4, 4] = 1 / G_TT
    S[5, 5] = 1 / G_L
    S[6, 6] = 1 / G_L
    return SMatrix{6,6}(inv(S))
end

# ──────────────────────────────────────────────────────────────────────────────
# Tensor ↔ Mandel / Voigt conversions
# ──────────────────────────────────────────────────────────────────────────────

"""
    voigt_to_4th(C6) -> Array{Float64,4}

Convert a 6×6 Voigt stiffness matrix to a full 3×3×3×3 fourth-order tensor.
"""
function voigt_to_4th(C6::AbstractMatrix)
    @assert size(C6) == (6, 6) "C6 must be a 6x6 matrix!"
    C4 = zeros(3, 3, 3, 3)
    vmap = [(1,1),(2,2),(3,3),(2,3),(1,3),(1,2)]
    for (α, (i,j)) in enumerate(vmap), (β, (k,l)) in enumerate(vmap)
        v = C6[α, β]
        C4[i,j,k,l] = v
        C4[j,i,k,l] = v
        C4[i,j,l,k] = v
        C4[j,i,l,k] = v
    end
    return C4
end


"""
    mandel_to_4th(C6m) -> Array{Float64,4}

Convert a 6×6 Mandel stiffness matrix to a 3×3×3×3 tensor.
"""
function mandel_to_4th(C6m::AbstractMatrix)
    @assert size(C6) == (6, 6) "C6 must be a 6x6 matrix!"
    vmap = [(1,1),(2,2),(3,3),(2,3),(1,3),(1,2)]
    C4 = zeros(3, 3, 3, 3)
    for (α, (i,j)) in enumerate(vmap), (β, (k,l)) in enumerate(vmap)
        fac = mandel_scale(i,j) * mandel_scale(k,l)
        v = C6m[α, β] / fac
        C4[i,j,k,l] = v
        C4[j,i,k,l] = v
        C4[i,j,l,k] = v
        C4[j,i,l,k] = v
    end
    return C4
end


"""
    fourth_to_voigt(C4) -> Matrix{Float64}

Convert a 3×3×3×3 fourth-order tensor to a 6×6 Voigt stiffness matrix.
"""
function fourth_to_voigt(C4::AbstractArray{<:Real,4})
    @assert size(C4) == (3, 3, 3, 3) "C4 must be size 3x3x3x3 !"
    vmap = [(1,1),(2,2),(3,3),(2,3),(1,3),(1,2)]
    C6 = zeros(6, 6)
    for (α, (i,j)) in enumerate(vmap), (β, (k,l)) in enumerate(vmap)
        C6[α, β] = C4[i,j,k,l]
    end
    return C6
end

"""
    fourth_to_mandel(C4) -> Matrix{Float64}

Convert a 3×3×3×3 fourth-order tensor to a 6×6 Mandel stiffness matrix.
"""
function fourth_to_mandel(C4::AbstractArray{<:Real,4})
    @assert size(C4) == (3, 3, 3, 3) "C4 must be size 3x3x3x3 !"
    vmap = [(1,1),(2,2),(3,3),(2,3),(1,3),(1,2)]
    C6 = zeros(6, 6)
    for (α, (i,j)) in enumerate(vmap), (β, (k,l)) in enumerate(vmap)
        fac = mandel_scale(i,j) * mandel_scale(k,l)
        C6[α, β] = fac * C4[i,j,k,l]
    end
    return C6
end



# ──────────────────────────────────────────────────────────────────────────────
# Fourth-order tensor algebra helpers
# ──────────────────────────────────────────────────────────────────────────────

"""Double contraction A:B of two 3×3×3×3 tensors → 3×3×3×3 tensor."""
function double_contract(A::AbstractArray{<:Real,4}, B::AbstractArray{<:Real,4})
    C = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = 0.0
        for m in 1:3, n in 1:3
            s += A[i,j,m,n] * B[m,n,k,l]
        end
        C[i,j,k,l] = s
    end
    return C
end

"""Fourth-order identity tensor I₄ : Iᵢⱼₖₗ = ½(δᵢₖδⱼₗ + δᵢₗδⱼₖ)"""
function identity4()
    I4 = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3
        I4[i,j,i,j] += 0.5
        I4[i,j,j,i] += 0.5
    end
    return I4
end

"""Invert a 3×3×3×3 tensor via its 6×6 Voigt representation."""
function invert4(A::AbstractArray{<:Real,4})
    A6 = fourth_to_voigt(A)
    return voigt_to_4th(inv(A6))
end

"""Add two 3×3×3×3 tensors with scalar prefactors: α*A + β*B"""
function add4(A, B; α=1.0, β=1.0)
    return α .* A .+ β .* B
end

# ──────────────────────────────────────────────────────────────────────────────
# Rotation of fourth-order stiffness tensor
# ──────────────────────────────────────────────────────────────────────────────

"""
    rotation_matrix(θ, ϕ=0.0, ψ=0.0) -> SMatrix{3,3}

Build a 3×3 rotation matrix from Euler angles (in radians).
Convention: Rz(ψ) · Ry(ϕ) · Rx(θ)  (extrinsic x-y-z).
For in-plane rotation about z-axis, use only θ.
"""
function rotation_matrix(θ::Real, ϕ::Real=0.0, ψ::Real=0.0)
    Rx = @SMatrix [1      0       0;
                   0  cos(θ)  -sin(θ);
                   0  sin(θ)   cos(θ)]
    Ry = @SMatrix [ cos(ϕ)  0  sin(ϕ);
                        0   1      0;
                   -sin(ϕ)  0  cos(ϕ)]
    Rz = @SMatrix [cos(ψ)  -sin(ψ)  0;
                   sin(ψ)   cos(ψ)  0;
                       0        0   1]
    return Rz * Ry * Rx
end

"""
    rotate_stiffness_4th(C4, Q) -> Array{Float64,4}

Rotate a 3×3×3×3 stiffness tensor by rotation matrix Q.
C'ᵢⱼₖₗ = Qᵢₘ Qⱼₙ Qₖₚ Qₗq Cₘₙₚq
"""
function rotate_stiffness_4th(C4::AbstractArray{<:Real,4}, Q::AbstractMatrix)
    C4r = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = 0.0
        for m in 1:3, n in 1:3, p in 1:3, q in 1:3
            s += Q[i,m] * Q[j,n] * Q[k,p] * Q[l,q] * C4[m,n,p,q]
        end
        C4r[i,j,k,l] = s
    end
    return C4r
end

"""
    rotate_stiffness_voigt(C6, θ, ϕ=0, ψ=0) -> Matrix{Float64}

Rotate a 6×6 Voigt stiffness matrix by Euler angles (radians).
"""
function rotate_stiffness_voigt(C6::AbstractMatrix, θ::Real, ϕ::Real=0.0, ψ::Real=0.0)
    Q = rotation_matrix(θ, ϕ, ψ)
    C4 = voigt_to_4th(C6)
    C4r = rotate_stiffness_4th(C4, Q)
    return fourth_to_voigt(C4r)
end

# ──────────────────────────────────────────────────────────────────────────────
# Eshelby tensors
# ──────────────────────────────────────────────────────────────────────────────

"""
    eshelby_tensor_needle(νm) -> Array{Float64,4}

Eshelby tensor for a needle-shaped (UD fiber, aspect ratio → ∞) inclusion
embedded in an isotropic matrix with Poisson ratio νm.

Components in Voigt notation (fiber along x₁):
Reference: Eshelby (1957), Tucker & Liang (1999).
"""
function eshelby_tensor_needle(νm::Real)
    ν = νm
    S = zeros(6, 6)
    # Non-zero Voigt components for needle (α → ∞)
    S[2,2] = S[3,3] = (5 - 4ν) / (8(1 - ν))
    S[2,3] = S[3,2] = (4ν - 1) / (8(1 - ν))
    S[4,4] = (3 - 4ν) / (4(1 - ν))   # S₂₃₂₃
    S[5,5] = S[6,6] = 1/4             # S₁₂₁₂ = S₁₃₁₃
    return voigt_to_4th(S)
end

"""
    eshelby_tensor_sphere(νm) -> Array{Float64,4}

Eshelby tensor for a spherical inclusion (aspect ratio = 1).
"""
function eshelby_tensor_sphere(νm::Real)
    ν = νm
    a = (7 - 5ν) / (15(1 - ν))
    b = (5ν - 1) / (15(1 - ν))
    c = (4 - 5ν) / (15(1 - ν))
    S = zeros(6, 6)
    S[1,1] = S[2,2] = S[3,3] = a
    S[1,2] = S[2,1] = S[1,3] = S[3,1] = S[2,3] = S[3,2] = b
    S[4,4] = S[5,5] = S[6,6] = c
    return voigt_to_4th(S)
end

"""
    eshelby_tensor_spheroid(νm, ar) -> Array{Float64,4}

Eshelby tensor for a prolate spheroidal inclusion with aspect ratio ar = a₁/a₂ > 1
(a₁ along x₁ = fiber axis). Uses closed-form expressions for prolate spheroids.

Reference: Mura (1987) "Micromechanics of Defects in Solids".
"""
function eshelby_tensor_spheroid(νm::Real, ar::Real)
    ν = νm
    α = ar  # a₁/a₂, with a₁ being the long axis

    if α ≈ 1.0
        return eshelby_tensor_sphere(ν)
    end

    if α > 1.0  # prolate
        g = α / (α^2 - 1)^(3/2) * (α * sqrt(α^2 - 1) - acosh(α))
    else  # oblate
        g = α / (1 - α^2)^(3/2) * (acos(α) - α * sqrt(1 - α^2))
    end

    S = zeros(6, 6)
    # Axial (11) components
    S[1,1] = 1 - (4(1-ν) + 1) / (2(1-ν)) * (1 - (2α^2*g)/(α^2-1)) / (2*(α^2-1))
    # These are the closed-form expressions for prolate spheroid
    # Using notation from Tucker & Liang (1999) / Advani & Tucker (1987)
    q = 2 / (α^2 - 1)

    I₁ = -g * 4π / (α^2 - 1) + 4π / (1 - 1/α^2)^(1/2) / α  # placeholder
    # Full closed-form Eshelby tensor for prolate spheroid:
    S11 = 1 / (2*(1-ν)) * (2*(1-ν) - (1+2α^2)/(α^2-1) * (1 - 2*α^2/(α^2-1)*(1-g)))
    S22_33 = 3*α^2 / (8*(1-ν)*(α^2-1)) + 1/(4*(1-ν)) * (1 - (1+3/2*α^2/(α^2-1))*g)
    S12_13 = -α^2 / (2*(1-ν)*(α^2-1)) + 1/(4*(1-ν)) * g*(1-1/(α^2-1)) * 3/2
    S21_31 = -1/(4*(1-ν)) + (α^2+1)/(4*(1-ν)*(α^2-1)) - g*(3*α^2-1)/(4*(1-ν)*(α^2-1))
    S23_32 = α^2/(4*(1-ν)*(α^2-1)) - 1/(4*(1-ν))*(g*(3*α^2-1)/(2*(α^2-1)) - 1)

    S44 = α^2/(2*(α^2-1)) + g*(1+α^2) / (4*(α^2-1)) - 1/4  # S₂₃₂₃ with ν correction
    S55_66 = 1/(4*(1-ν))*((α^2/(2*(α^2-1)) + 1 - ν) - g*(α^2/(α^2-1)/2 + 1 - 2ν))

    # Re-derive carefully using Mura's book notation
    # a = a₁ (long), b = a₂ = a₃ (short); α = a/b
    a2 = α^2
    if α > 1.0
        Iₐ = 4π*α / (a2-1)^(3/2) * (α/sqrt(a2-1) - acosh(α))
    else
        Iₐ = 4π*α / (1-a2)^(3/2) * (acos(α) - α*sqrt(1-a2))
    end
    Ib = (4π - Iₐ) / 2

    Iₐₐ = (4π/3 - Iₐ) / (a2-1)
    Iab = (Ib - Iₐ) / (a2-1)
    Ibb = (4π/3 - 2*Iab - Iₐₐ) / 3   # Only one distinct transverse component needed

    # Eshelby components (Mura 1987 Table 11.1)
    S[1,1] = 3/(8π*(1-ν)) * a2*Iₐₐ + (1-2ν)/(4π*(1-ν)) * Iₐ
    S[2,2] = 3/(8π*(1-ν)) * Ibb + (1-2ν)/(4π*(1-ν)) * Ib
    S[3,3] = S[2,2]
    S[1,2] = 1/(8π*(1-ν)) * Iab - (1-2ν)/(8π*(1-ν)) * Iₐ
    S[1,3] = S[1,2]
    S[2,1] = 1/(8π*(1-ν)) * a2*Iₐₐ - (1-2ν)/(8π*(1-ν)) * Ib  # S₂₂₁₁ type swap
    S[3,1] = S[2,1]
    S[2,3] = 1/(8π*(1-ν)) * Ibb - (1-2ν)/(8π*(1-ν)) * Ib
    S[3,2] = S[2,3]
    S[4,4] = (a2+1)/(16π*(1-ν)) * (2*Iab + 2*Ibb) / 2 + (1-2ν)/(8π*(1-ν)) * (Ib + Iₐ) / 2  # approx
    # Shear components
    Sab = (a2+1)/(16π*(1-ν)) * Iab - (1-2ν)/(8π*(1-ν)) * (Iₐ+Ib)/4  # S₁₂₁₂ type
    S[5,5] = Sab
    S[6,6] = Sab
    S44v = 1/(16π*(1-ν))*(Ibb+Ibb)/2*(2) + (1-2ν)/(8π*(1-ν))*(Ib+Ib)/4  # S₂₃₂₃
    S[4,4] = 1/(16π*(1-ν)) * 2*Ibb + (1-2ν)/(8π*(1-ν)) * 2*Ib/2

    # Corrections: use simpler Tucker & Liang parametric form
    return _eshelby_prolate_TL(νm, α)
end

"""Internal: Eshelby tensor for prolate spheroid via Tucker & Liang (1999) expressions."""
function _eshelby_prolate_TL(ν::Real, α::Real)
    if α ≈ 1.0
        return eshelby_tensor_sphere(ν)
    end
    a = α
    if a > 1.0
        g = a/(a^2-1)^(3/2) * (a*sqrt(a^2-1) - acosh(a))
    else
        g = a/(1-a^2)^(3/2) * (acos(a) - a*sqrt(1-a^2))
    end

    S = zeros(6, 6)
    S[1,1] = 1 - 2ν + (3a^2-1)/(a^2-1) - g*(1-2ν + 3a^2/(a^2-1))
    S[1,1] /= (2*(1-ν))

    S[2,2] = 3*a^2/(8*(1-ν)*(a^2-1)) + (1-2ν)/(4*(1-ν)) * (1 - (1 + 3a^2/(2*(a^2-1)))*g)
    S[3,3] = S[2,2]

    S[1,2] = -a^2/(2*(1-ν)*(a^2-1)) + (1/(4*(1-ν))) * g * (3a^2/(2*(a^2-1)) - (1-2ν))
    S[1,3] = S[1,2]

    S[2,1] = -1/(2*(1-ν)) * (1 - (a^2+1)/(a^2-1) + g*(1-2ν + 3/(2*(a^2-1))))
    # correct S21: should be 1/(2*(1-ν)) * (-1 + a^2/(a^2-1)... )
    S[2,1] = 1/(2*(1-ν)) * (-(1-2ν) + 1/(2*(a^2-1)) - g*(3/(2*(a^2-1)) - (1-2ν)))
    S[3,1] = S[2,1]

    S[2,3] = a^2/(4*(1-ν)*(a^2-1)) - 1/(4*(1-ν)) * (1 + (3*a^2-1)/(a^2-1)*g/2 - 1)
    S[3,2] = S[2,3]

    # Shear components
    S[4,4] = (a^2/(2*(a^2-1)) + (1-2ν)/4 - g*(3*a^2/(4*(a^2-1)) + (1-2ν)/4)) / (1-ν)
    # Corrected: S₂₃₂₃
    S[4,4] = 1/(4*(1-ν)) * (a^2/(a^2-1)/2 + (1-2ν) - g*(3a^2/(2*(a^2-1)) - (1-2ν)))
    # S₁₂₁₂ = S₁₃₁₃
    S[5,5] = 1/(4*(1-ν)) * (-(1-2ν) + 1/(2*(a^2-1)) + g*((1+a^2)/(2*(a^2-1)) + (1-2ν)))
    S[6,6] = S[5,5]

    return voigt_to_4th(S)
end

# ──────────────────────────────────────────────────────────────────────────────
# Orientation averaging (Advani-Tucker)
# ──────────────────────────────────────────────────────────────────────────────

"""
    orientation_average(C4_UD, A2, A4) -> Array{Float64,4}

Orientation-average a UD stiffness tensor C4_UD using second-order fiber
orientation tensor A2 (3×3) and fourth-order orientation tensor A4 (3×3×3×3).

Uses the Advani-Tucker (1987) closure approximation:
Cᵢⱼₖₗ = B₁ A₄ᵢⱼₖₗ + B₂(A₂ᵢⱼδₖₗ + A₂ₖₗδᵢⱼ) + B₃(A₂ᵢₖδⱼₗ + A₂ᵢₗδⱼₖ + A₂ⱼₖδᵢₗ + A₂ⱼₗδᵢₖ) + B₄δᵢⱼδₖₗ + B₅(δᵢₖδⱼₗ + δᵢₗδⱼₖ)

where B coefficients are derived from the UD stiffness invariants.
"""
function orientation_average(C4_UD::AbstractArray{<:Real,4}, A2::AbstractMatrix, A4::AbstractArray{<:Real,4})
    # Extract independent UD stiffness components (Voigt, fiber along x₁)
    C6 = fourth_to_voigt(C4_UD)
    c11 = C6[1,1]; c22 = C6[2,2]; c12 = C6[1,2]; c23 = C6[2,3]; c44 = C6[4,4]; c55 = C6[5,5]

    # Advani-Tucker B-coefficients
    B1 = c11 + c22 - 2*c12 - 4*c55
    B2 = c12 - c23
    B3 = c55
    B4 = c23
    B5 = (c22 - c23) / 2

    δ = I(3)
    C4 = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        val = 0.0
        val += B1 * A4[i,j,k,l]
        val += B2 * (A2[i,j]*δ[k,l] + A2[k,l]*δ[i,j])
        val += B3 * (A2[i,k]*δ[j,l] + A2[i,l]*δ[j,k] + A2[j,k]*δ[i,l] + A2[j,l]*δ[i,k])
        val += B4 * δ[i,j]*δ[k,l]
        val += B5 * (δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k])
        C4[i,j,k,l] = val
    end
    return C4
end

"""
    a4_linear_closure(A2) -> Array{Float64,4}

Linear closure approximation for the fourth-order orientation tensor from A2.
"""
function a4_linear_closure(A2::AbstractMatrix)
    δ = I(3)
    A4 = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        A4[i,j,k,l] = (-1/24) * (δ[i,j]*δ[k,l] + δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k]) +
                       (1/6)  * (A2[i,j]*δ[k,l] + A2[i,k]*δ[j,l] + A2[i,l]*δ[j,k] +
                                  A2[k,l]*δ[i,j] + A2[j,l]*δ[i,k] + A2[j,k]*δ[i,l])
    end
    return A4
end

"""
    a4_quadratic_closure(A2) -> Array{Float64,4}

Quadratic (IBOF) closure approximation for the fourth-order orientation tensor.
"""
a4_quadratic_closure(A2::AbstractMatrix) = [A2[i,j]*A2[k,l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

# ──────────────────────────────────────────────────────────────────────────────
# Mori-Tanaka scheme
# ──────────────────────────────────────────────────────────────────────────────

"""
    dilute_strain_concentration(Cf, Cm, S4) -> Array{Float64,4}

Compute the dilute strain concentration (A-tensor) for inclusions.
A = [I + S : C_m⁻¹ : (C_f - C_m)]⁻¹
"""
function dilute_strain_concentration(Cf::AbstractArray{<:Real,4},
                                      Cm::AbstractArray{<:Real,4},
                                      S4::AbstractArray{<:Real,4})
    I4 = identity4()
    Cm_inv = invert4(Cm)
    diff = add4(Cf, Cm; α=1.0, β=-1.0)
    inner = add4(I4, double_contract(S4, double_contract(Cm_inv, diff)))
    return invert4(inner)
end

"""
    mori_tanaka(Cm6, Cf6, vf, eshelby_S4; symmetrize=false) -> Matrix{Float64}

Mori-Tanaka effective stiffness (Benveniste formulation) for a two-phase
composite.

Arguments:
  Cm6       — 6×6 Voigt stiffness of the matrix
  Cf6       — 6×6 Voigt stiffness of the filler/fiber (UD orientation: x₁)
  vf        — fiber volume fraction (0–1)
  eshelby_S4 — 3×3×3×3 Eshelby tensor for the inclusion shape
  symmetrize — if true, apply Segura (2023) symmetrization

Returns the 6×6 effective Voigt stiffness tensor.
"""
function mori_tanaka(Cm6::AbstractMatrix, Cf6::AbstractMatrix, vf::Real,
                     eshelby_S4::AbstractArray{<:Real,4};
                     symmetrize::Bool=false)
    vm = 1 - vf
    Cm4 = voigt_to_4th(Cm6)
    Cf4 = voigt_to_4th(Cf6)
    I4  = identity4()

    A4  = dilute_strain_concentration(Cf4, Cm4, eshelby_S4)

    # Mori-Tanaka concentration tensor B = A · [vm·I + vf·A]⁻¹
    inner = add4(I4, A4; α=vm, β=vf)
    B4    = double_contract(A4, invert4(inner))

    # Effective stiffness: C* = Cm + vf·(Cf - Cm)·B
    diff  = add4(Cf4, Cm4; α=1.0, β=-1.0)
    C_eff = add4(Cm4, double_contract(diff, B4); α=1.0, β=vf)

    if symmetrize
        C_eff = _symmetrize4(C_eff)
    end
    return fourth_to_voigt(C_eff)
end

"""
    _symmetrize4(C4) -> Array{Float64,4}

Enforce major symmetry Cᵢⱼₖₗ = Cₖₗᵢⱼ via simple averaging.
(Segura et al. 2023 correction for non-symmetric Mori-Tanaka results.)
"""
function _symmetrize4(C4::AbstractArray{<:Real,4})
    C4s = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        C4s[i,j,k,l] = 0.5 * (C4[i,j,k,l] + C4[k,l,i,j])
    end
    return C4s
end

"""
    mori_tanaka_multi(Cm6, phases; symmetrize=false) -> Matrix{Float64}

Multi-phase Mori-Tanaka for N inclusion types.

`phases` is a vector of named tuples:
  (C6=<6×6 stiffness>, vf=<volume fraction>, S4=<Eshelby tensor>)

The matrix volume fraction is computed as 1 - Σvf_i.
"""
function mori_tanaka_multi(Cm6::AbstractMatrix, phases::AbstractVector;
                           symmetrize::Bool=false)
    Cm4 = voigt_to_4th(Cm6)
    I4  = identity4()
    vm  = 1.0 - sum(p.vf for p in phases)

    # Numerator and denominator of the MT formula
    num = deepcopy(Cm4) .* vm
    den = I4 .* vm

    for p in phases
        Cf4 = voigt_to_4th(p.C6)
        A4  = dilute_strain_concentration(Cf4, Cm4, p.S4)
        num = add4(num, double_contract(Cf4, A4); α=1.0, β=p.vf)
        den = add4(den, A4; α=1.0, β=p.vf)
    end

    C_eff = double_contract(num, invert4(den))

    if symmetrize
        C_eff = _symmetrize4(C_eff)
    end
    return fourth_to_voigt(C_eff)
end

"""
    mori_tanaka_orientation(Cm6, Cf6, vf, eshelby_S4, A2;
                            closure=:linear, symmetrize=false) -> Matrix{Float64}

Mori-Tanaka with orientation averaging (Advani-Tucker 1987).

Computes the UD (aligned) effective stiffness, then applies orientation
averaging using second-order fiber orientation tensor A2.

Arguments:
  A2      — 3×3 fiber orientation tensor (must satisfy tr(A2) = 1)
  closure — :linear, :quadratic (closure approx for A4)
"""
function mori_tanaka_orientation(Cm6::AbstractMatrix, Cf6::AbstractMatrix,
                                 vf::Real, eshelby_S4::AbstractArray{<:Real,4},
                                 A2::AbstractMatrix;
                                 closure::Symbol=:linear,
                                 symmetrize::Bool=false)
    # Step 1: UD effective stiffness (fiber along x₁)
    C_UD6 = mori_tanaka(Cm6, Cf6, vf, eshelby_S4; symmetrize=symmetrize)
    C_UD4 = voigt_to_4th(C_UD6)

    # Step 2: Fourth-order orientation tensor from A2
    A4 = if closure == :quadratic
        a4_quadratic_closure(A2)
    else
        a4_linear_closure(A2)
    end

    # Step 3: Orientation average
    C_avg4 = orientation_average(C_UD4, A2, A4)
    return fourth_to_voigt(C_avg4)
end

# ──────────────────────────────────────────────────────────────────────────────
# Halpin-Tsai model (scalar, with optional Cox Shear-Lag)
# ──────────────────────────────────────────────────────────────────────────────

"""
    halpin_tsai_param(Ef, Em, ξ) -> (η, modulus)

Core Halpin-Tsai interpolation parameter η for a given reinforcing modulus Ef,
matrix modulus Em, and shape/geometry parameter ξ.

Returns (η, composite modulus E).
"""
function halpin_tsai_param(Ef::Real, Em::Real, vf::Real, ξ::Real)
    η = (Ef/Em - 1) / (Ef/Em + ξ)
    E = Em * (1 + ξ * η * vf) / (1 - η * vf)
    return η, E
end

"""
    halpin_tsai_UD(Em, νm, Ef, νf, Gf, vf, ar; use_shear_lag=false, l_fiber=nothing, d_fiber=nothing) -> Matrix{Float64}

Halpin-Tsai homogenization for a UD fiber-reinforced polymer.

Returns the 6×6 Voigt stiffness matrix of the transversely isotropic composite.

Arguments:
  Em, νm   — matrix Young's modulus and Poisson ratio
  Ef, νf   — fiber Young's modulus and Poisson ratio (longitudinal / transverse for TI fiber)
  Gf       — fiber shear modulus
  vf       — fiber volume fraction
  ar       — fiber aspect ratio (length / diameter)
  use_shear_lag — apply Cox (1952) Shear-Lag correction to longitudinal modulus
  l_fiber  — fiber length [m] (required if use_shear_lag=true)
  d_fiber  — fiber diameter [m] (required if use_shear_lag=true)
"""
function halpin_tsai_UD(Em::Real, νm::Real, Ef::Real, νf::Real, Gf::Real, vf::Real, ar::Real;
                        use_shear_lag::Bool=false, l_fiber=nothing, d_fiber=nothing)
    Gm = Em / (2*(1 + νm))
    vm = 1 - vf

    # Longitudinal modulus (rule of mixtures)
    E_L = Ef * vf + Em * vm

    # Shear-Lag correction (Cox 1952)
    if use_shear_lag
        if isnothing(l_fiber) || isnothing(d_fiber)
            error("l_fiber and d_fiber must be provided when use_shear_lag=true")
        end
        # Cox factor β
        β = sqrt(2π * Gm / (Ef * log(1/vf))) / (d_fiber / 2)
        βL = β * l_fiber / 2
        η_cox = 1 - tanh(βL) / βL
        E_L = Ef * η_cox * vf + Em * vm
    end

    # Transverse modulus (Halpin-Tsai, ξ = 2 for circular fiber)
    ξT = 2.0
    _, E_T = halpin_tsai_param(Ef, Em, vf, ξT)

    # In-plane shear modulus (Halpin-Tsai, ξ = 1)
    ξG = 1.0
    _, G_L = halpin_tsai_param(Gf, Gm, vf, ξG)

    # Transverse shear modulus (Halpin-Tsai)
    Gf_T = Gf  # isotropic fiber assumption
    _, G_TT = halpin_tsai_param(Gf_T, Gm, vf, ξG)

    # Longitudinal Poisson ratio (rule of mixtures)
    ν_L = νf * vf + νm * vm

    # Transverse Poisson ratio from transverse isotropy constraint
    ν_TT = E_T / (2 * G_TT) - 1

    return transverse_isotropic_stiffness(E_L, E_T, ν_L, ν_TT, G_L)
end

"""
    halpin_tsai(Em, νm, Ef, νf, Gf, vf, ar, θ_distribution;
                use_shear_lag=false, l_fiber=nothing, d_fiber=nothing) -> Matrix{Float64}

Halpin-Tsai with in-plane fiber angle distribution via laminate analogy.

`θ_distribution` is a vector of (angle_rad, weight) tuples representing the
fiber orientation distribution function (ODF) discretization.

The composite stiffness is the weighted average of rotated UD stiffnesses.
"""
function halpin_tsai(Em::Real, νm::Real, Ef::Real, νf::Real, Gf::Real, vf::Real, ar::Real,
                     θ_distribution::AbstractVector{<:Tuple};
                     use_shear_lag::Bool=false, l_fiber=nothing, d_fiber=nothing)
    C_UD6 = halpin_tsai_UD(Em, νm, Ef, νf, Gf, vf, ar;
                           use_shear_lag=use_shear_lag, l_fiber=l_fiber, d_fiber=d_fiber)

    C_eff = zeros(6, 6)
    total_weight = sum(w for (_, w) in θ_distribution)

    for (θ, w) in θ_distribution
        C_rot = rotate_stiffness_voigt(C_UD6, 0.0, 0.0, θ)  # rotate about z-axis
        C_eff .+= (w / total_weight) .* C_rot
    end
    return C_eff
end

# ──────────────────────────────────────────────────────────────────────────────
# Material property extraction
# ──────────────────────────────────────────────────────────────────────────────

"""
    compliance_from_stiffness(C6) -> Matrix{Float64}

Compute the 6×6 compliance matrix S = C⁻¹.
"""
compliance_from_stiffness(C6::AbstractMatrix) = inv(C6)

"""
    young_modulus_from_stiffness(C6, direction) -> Float64

Extract the directional Young's modulus from a 6×6 stiffness matrix.
`direction` is a unit vector (3-element). Uses S = C⁻¹ and projects.
"""
function young_modulus_from_stiffness(C6::AbstractMatrix, direction::AbstractVector)
    S6 = compliance_from_stiffness(C6)
    n = normalize(direction)
    # Young's modulus: 1/E = n⊗n : S : n⊗n  (contracted)
    # In Voigt: 1/E = Σ nᵢnⱼnₖnₗ Sᵢⱼₖₗ
    S4 = voigt_to_4th(S6)
    inv_E = 0.0
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        inv_E += n[i]*n[j]*n[k]*n[l] * S4[i,j,k,l]
    end
    return 1 / inv_E
end

"""
    young_moduli(C6) -> (E1, E2, E3)

Extract principal Young's moduli E₁, E₂, E₃ from compliance S = C⁻¹.
"""
function young_moduli(C6::AbstractMatrix)
    S6 = compliance_from_stiffness(C6)
    return 1/S6[1,1], 1/S6[2,2], 1/S6[3,3]
end

"""
    shear_moduli(C6) -> (G23, G13, G12)

Extract shear moduli from compliance S = C⁻¹.
"""
function shear_moduli(C6::AbstractMatrix)
    S6 = compliance_from_stiffness(C6)
    return 1/S6[4,4], 1/S6[5,5], 1/S6[6,6]
end

"""
    poisson_ratios(C6) -> (ν12, ν13, ν23)

Extract major Poisson ratios from compliance.
"""
function poisson_ratios(C6::AbstractMatrix)
    S6 = compliance_from_stiffness(C6)
    ν12 = -S6[1,2] * (1/S6[1,1])
    ν13 = -S6[1,3] * (1/S6[1,1])
    ν23 = -S6[2,3] * (1/S6[2,2])
    return ν12, ν13, ν23
end

# ──────────────────────────────────────────────────────────────────────────────
# Orientation tensor utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
    a2_isotropic() -> Matrix{Float64}

Second-order orientation tensor for a fully isotropic (3D random) distribution.
A2 = I/3
"""
a2_isotropic() = Matrix(Diagonal([1/3, 1/3, 1/3]))

"""
    a2_planar_isotropic() -> Matrix{Float64}

Second-order orientation tensor for a 2D planar isotropic distribution
(random in x₁-x₂ plane, fibers lie in plane).
"""
a2_planar_isotropic() = Matrix(Diagonal([1/2, 1/2, 0.0]))

"""
    a2_unidirectional(; axis=1) -> Matrix{Float64}

Second-order orientation tensor for perfectly aligned fibers along `axis` (1, 2, or 3).
"""
function a2_unidirectional(; axis::Int=1)
    A2 = zeros(3, 3)
    A2[axis, axis] = 1.0
    return A2
end

"""
    a2_from_angles(θs, ϕs, weights) -> Matrix{Float64}

Build A2 from a discrete set of fiber orientations (θ polar, ϕ azimuthal) with weights.
"""
function a2_from_angles(θs::AbstractVector, ϕs::AbstractVector, weights::AbstractVector)
    @assert length(θs) == length(ϕs) == length(weights)
    A2 = zeros(3, 3)
    W  = sum(weights)
    for (θ, ϕ, w) in zip(θs, ϕs, weights)
        p = @SVector [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
        for i in 1:3, j in 1:3
            A2[i,j] += (w/W) * p[i] * p[j]
        end
    end
    return A2
end

end # module HomoPy
using .HomoPy
# ──────────────────────────────────────────────────────────────────────────────
# Example / demonstration (run as a script)
# ──────────────────────────────────────────────────────────────────────────────
function run_me()
   

    println("="^60)
    println("HomoPy.jl — Demonstration")
    println("="^60)

    # Material parameters (typical CF/PA6 values)
    Em, νm = 3.35, 0.385          # PA6 matrix
    Ef, νf = 230, 0.26            # Carbon fiber (longitudinal)
    Gf     = Ef / (2*(1 + νf))
    vf     = 0.20                   # 20% fiber volume fraction
    ar     = 20.0                   # aspect ratio

    println("\n--- Halpin-Tsai (UD, no Shear-Lag) ---")
    C_HT = HomoPy.halpin_tsai_UD(Em, νm, Ef, νf, Gf, vf, ar)
    E1, E2, E3 = HomoPy.young_moduli(C_HT)
    G23, G13, G12 = HomoPy.shear_moduli(C_HT)
    ν12, ν13, ν23 = HomoPy.poisson_ratios(C_HT)
    println("  E₁ = $(round(E1 , digits=2)) GPa")
    println("  E₂ = $(round(E2 , digits=2)) GPa")
    println("  G₁₂ = $(round(G12 , digits=3)) GPa")
    println("  ν₁₂ = $(round(ν12, digits=4))")

    println("\n--- Halpin-Tsai (UD, with Cox Shear-Lag) ---")
    C_HT_SL = HomoPy.halpin_tsai_UD(Em, νm, Ef, νf, Gf, vf, ar;
                                     use_shear_lag=true, l_fiber=1e-3, d_fiber=7e-6)
    E1_SL, _, _ = HomoPy.young_moduli(C_HT_SL)
    println("  E₁ (Shear-Lag) = $(round(E1_SL , digits=2)) GPa")

    println("\n--- Mori-Tanaka (UD needle inclusions) ---")
    Cm6 = HomoPy.isotropic_stiffness(Em, νm)
    Cf6 = HomoPy.isotropic_stiffness(Ef, νf)   # isotropic fiber approximation
    S4  = HomoPy.eshelby_tensor_needle(νm)
    C_MT6 = HomoPy.mori_tanaka(Cm6, Cf6, vf, S4; symmetrize=false)
    E1_MT, E2_MT, _ = HomoPy.young_moduli(C_MT6)
    println("  E₁ = $(round(E1_MT , digits=2)) GPa")
    println("  E₂ = $(round(E2_MT , digits=2)) GPa")

    println("\n--- Mori-Tanaka with orientation averaging (planar isotropic) ---")
    A2 = HomoPy.a2_planar_isotropic()
    C_MT_OA = HomoPy.mori_tanaka_orientation(Cm6, Cf6, vf, S4, A2; closure=:linear)
    E1_OA, E2_OA, E3_OA = HomoPy.young_moduli(C_MT_OA)
    println("  E₁ = $(round(E1_OA , digits=2)) GPa")
    println("  E₂ = $(round(E2_OA , digits=2)) GPa")
    println("  E₃ = $(round(E3_OA , digits=2)) GPa")

    println("\n--- Symmetry check (major symmetry Cᵢⱼₖₗ = Cₖₗᵢⱼ) ---")
    C4_MT = HomoPy.voigt_to_4th(C_MT6)
    max_asym = maximum(abs(C4_MT[i,j,k,l] - C4_MT[k,l,i,j])
                       for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    println("  max|Cᵢⱼₖₗ - Cₖₗᵢⱼ| = $(round(max_asym, sigdigits=3)) Pa")


    println("\nMulti Phase\n")

    #matrix
    E_m , nu_m = 2.0, 0.4
    C_m = HomoPy.isotropic_stiffness(E_m, nu_m)

    #carbon vf = 0.22
    E1_carbon, E2_carbon, nu12_carbon, nu23_carbon, G12_carbon = 230.0, 20.0, 0.07, 0.39, 40.0
    C_carbon = HomoPy.transverse_isotropic_stiffness(E1_carbon, E2_carbon, nu12_carbon, nu23_carbon, G12_carbon )
    ar_carbon = 22.0
    vf_carbon = 0.229
    S_carbon = eshelby_tensor_spheroid(nu_m, ar_carbon)


    #glass
    E_glass, nu_glass = 75.0, 0.23
    C_glass = isotropic_stiffness(E_glass, nu_glass)
    vf_glass = 0.0826
    ar_glass = 25.0
    S_glass = eshelby_tensor_spheroid(nu_m, ar_glass)
    #(C6=<6×6 stiffness>, vf=<volume fraction>, S4=<Eshelby tensor>)
    phases = [
        (;C6 = C_carbon, vf = vf_carbon, S4 = S_carbon),
        (;C6 = C_glass, vf = vf_glass, S4 = S_glass),
    ]
    Ceff = mori_tanaka_multi(C_m, phases; symmetrize = true)               
    @show Ceff

    @show E1_eff, E2_eff, E3_eff = HomoPy.young_moduli(Ceff)


    println("\nDone.")
    return Ceff
end
