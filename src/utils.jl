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