#kron delta
@inline function δ(i, j) 
    i == j ? 1 : 0
end



"""
    convert_3333_to_66(tens)

    returns 4th-order tensor in engineering notation
"""
function convert_3333_to_66(tens; mandel = false)
    @assert size(tens) == (3,3,3,3)
    b = mandel ? sqrt(2) : 1
    
    code = ((1,1), (2,2), (3,3), (2,3), (1,3), (1, 2))
    m = (1, 1, 1, b, b, b)

    return SMatrix{6,6}(m[i] * m[j] * tens[code[i]..., code[j]...] for i in 1:6, j in 1:6)
end

function convert_66_to_3333(C66; mandel = false)

    lookup = (
        (1, 6, 5),
        (6, 2, 4),
        (5, 4, 3)
    )

    T = eltype(C66)
    # Helper to get the scaling factor based on notation
    function get_scale(α, mandel)
        !mandel && return one(T)
        if α > 3 
            return T(1/sqrt(2))
        end
        return one(T)
    end


    return SymmetricTensor{4,3}(
        (i, j, k, l) -> begin
                α = lookup[i][j]
                β = lookup[k][l]
                val = C66[α, β] * get_scale(α, mandel) * get_scale(β, mandel)
            end
    )

end






"""
    to_volume_fractions(weight_fractions::AbstractVector{T}, densities::AbstractVector{T}) where T<:Real

Convert weight fractions to volume fractions
Inputs:
    weight_fractions - AbstractVector of length 1 less than densities. Sum of weight_fractions must be between 0 and 1
    densities = AbstractVector of Reals

    Assumes that the 1st density is that of the matrix, whose weight fraction is (1-Σwi)

Ouputs:
    volume_fractions - Vector{T} of volume_fractions with length of densities. The 1st values is that of the matrix.
"""
function to_volume_fractions(weight_fractions::AbstractVector{T}, densities::AbstractVector{T}) where T<:Real

    @assert all(weight_fractions .> 0) "Weigth fractions must be positive!"
    @assert all(densities .> 0) "Densities must be positive!" 
    @assert length(weight_fractions) == length(densities) - 1 "Length of weight fractions vector must 1 less than densities vector!"

    N = length(densities)
    W = sum(weight_fractions)
    @assert 0 <= W <= 1 "Sum of weight fractions must be between 0 and 1!"

    D = (1- W) / densities[1] #start with the matrix contribution
    volume_fractions = zeros(T, length(densities))
    volume_fractions[1] =  D #matrix

    for (i,(w, rho)) in enumerate(zip(weight_fractions, view(densities, 2:N)))
        vi = w / rho
        volume_fractions[i + 1] = vi
        D += vi
    end
    volume_fractions ./= D
    return volume_fractions
end



"""
    to_weight_fractions(volume_fractions::AbstractVector{T}, densities::AbstractVector{T}) where T<:Real


Convert volume fractions to weight fractions.
Inputs:
    volume_fractions - AbstractVector of Reals, with length 1 less than densities.
    densities - AbstractVector of Reals, 1st value is assumed matrix.

Outputs:
    weight_fractions - Vector{T} with length of densities vector. First value is that of the matrix.
"""
function to_weight_fractions(volume_fractions::AbstractVector{T}, densities::AbstractVector{T}) where T<:Real
    @assert all(volume_fractions .> 0) "Volume fractions must be positive!"
    @assert all(densities .> 0) "Densities must be positive!" 
    @assert length(volume_fractions) == length(densities)-1 "Length of volume fraction vector must 1 less than densities vector!"

    N = length(densities)
    V = sum(volume_fractions)
    @assert 0<=V<=1 "Sum of volume fraction must be between 0 and 1!"

    D = densities[1] * (1 - V)
    weight_fractions = zeros(T, length(densities))
    weight_fractions[1] = D

    for (i, (v, rho)) in enumerate(zip(volume_fractions, view(densities, 2:N)))
        wi = rho * v
        D += wi
        weight_fractions[i + 1] = wi
    end
    weight_fractions ./= D
    return weight_fractions
end