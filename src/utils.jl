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

    # Helper to get the scaling factor based on notation
    function get_scale(α, mandel)
        if !mandel
            return 1.0 # Standard Voigt: scaling handled elsewhere or not needed depending on C vs S
        end
        # Mandel scaling: factors of sqrt(2) for indices 4, 5, 6
        s = 1.0
        if α > 3; s *= 1/sqrt(2); end
        return s
    end

    # Construct the 3x3x3x3 tensor
    # If mandel=true, we divide the C66 components by sqrt(2) for each shear index
    C4 = @SArray [
        begin
            α = lookup[i][j]
            β = lookup[k][l]
            val = C66[α, β]
            mandel ? val * get_scale(α, mandel) * get_scale(β, mandel) : val
        end
        for i=1:3, j=1:3, k=1:3, l=1:3
    ]
    
    return C4
end