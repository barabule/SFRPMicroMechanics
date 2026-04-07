#####ORIENTATION TENSORS 


abstract type AbstractOrientationTensor end


struct OrientationTensor{T<:Real} <:AbstractOrientationTensor
    a11::T
    a22::T

    function OrientationTensor(a11, a22)
        args = promote(a11, a22)
        T = eltype(args)
        a11, a22 = args[1], args[2]
        a11 = clamp(a11, T(1/3), T(1))
        a22 = clamp(a22, T(1/2 * (1 - a11)), min(a11, T(1-a11)))
        if a22 > a11 #just to make sure
            a11, a22 = a22, a11
        end
        return new{T}(a11, a22)
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
        #special case
        a1, a2, a4, a5, a6 = args    
        #a11 and a22 must be positive 
        
        a3 = 1 -a1 -a2
        m = SMatrix{3,3}([a1 a6 a5;
                          a6 a2 a4;
                          a5 a4 a3])
        λ, V = eigen(m)
        Σ = sum(λ)
        λ = λ * 1/Σ
        Λ = SMatrix{3,3}([λ[1] 0 0;
                          0 λ[2] 0;
                          0 0  λ[3]])
        m_new = V * Λ * V'
        a1, a2, a4, a5, a6 = m_new[1,1], m_new[2,2], m_new[2,3], m_new[1,3], m_new[1,2] 
       

        new{T}(a1, a2, a4, a5, a6)
    end
end

function get_all_coefficients(a::FullOrientationTensor)
    a3 = 1 - a.a1 - a.a2
    return (a.a1, a.a2, a3, a.a4, a.a5, a.a6)
end


function FullOrientationTensor(;a11 = nothing, a22 = nothing, a23 = nothing, a13 = nothing, a12 = nothing)
    @assert !(isnothing(a11) || isnothing(a22) || isnothing(a23) || isnothing(a13) || isnothing(a12))
    a1, a2, a4, a5, a6 = a11, a22, a23, a13, a12
    return FullOrientationTensor(a1, a2, a4, a5, a6)
end



function decompose(a::AbstractOrientationTensor)
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
        println(io, "A11 = $(p.a11)")
        println(io, "A22 = $(p.a22)")
        println(io, "A33 = $(1 - p.a11 - p.a22)")
    elseif isa(p, FullOrientationTensor)
        println(io, "A11 = $(p.a1)")
        println(io, "A22 = $(p.a2)")
        println(io, "A33 = $(1 - p.a1 - p.a2)")
        println(io, "A23 = $(p.a4)")
        println(io, "A13 = $(p.a5)")
        println(io, "A12 = $(p.a6)")
    end
end



function to_matrix(a::AbstractOrientationTensor)
    a1, a2, a3, a4, a5,a6  = get_all_coefficients(a)
    return SymmetricTensor{2, 3}(
                    [a1 a6 a5;
                     a6 a2 a4;
                     a5 a4 a3])
end