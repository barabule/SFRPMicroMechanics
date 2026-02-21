
"""
    tens2eng(tens)

    returns 4th-order tensor in engineering notation
"""
function tens2eng(tens; mandel = false)
    T = eltype(tens)
    A = MMatrix{6,6, T}(zeros(6, 6))
    b = mandel ? sqrt(2) : 1
    c = b
    # code = ((1,1), (2,2), (3,3), (2,3), (1,3), (1, 2))
    # for i in 1:6
    #     for j in 1:6
    #         A[i, j] = tens[code[i]..., code[j]...]
    #     end
    # end
    g =tens
    A = SMatrix{6,6,T}([g[1,1,1,1]     g[1,1,2,2]     g[1,1,3,3]     b * g[1,1,2,3]     b * g[1,1,1,3]     b * g[1,1,1,2];
                        g[2,2,1,1]     g[2,2,2,2]     g[2,2,3,3]     b * g[2,2,2,3]     b * g[2,2,1,3]     b * g[2,2,1,2];
                        g[3,3,1,1]     g[3,3,2,2]     g[3,3,3,3]     b * g[3,3,2,3]     b * g[3,3,1,3]     b * g[3,3,1,2];
                    c * g[2,3,1,1] c * g[2,3,2,2] c * g[2,3,3,3] b * c * g[2,3,2,3] b * c * g[2,3,1,3] b * c * g[2,3,1,2];
                    c * g[1,3,1,1] c * g[1,3,2,2] c * g[1,3,3,3] b * c * g[1,3,2,3] b * c * g[1,3,1,3] b * c * g[1,3,1,2];
                    c * g[1,2,1,1] c * g[1,2,2,2] c * g[1,2,3,3] b * c * g[1,2,2,3] b * c * g[1,2,1,3] b * c * g[1,2,1,2]] )
    return A
end

function convert_voigt_to_tensor(C66)

    lookup = (
        (1, 6, 5),
        (6, 2, 4),
        (5, 4, 3)
    )

    # Generate the 4th order tensor using a comprehension
    # This creates an SArray{Tuple{3,3,3,3}}
    C4 = @SArray [C66[lookup[i][j], lookup[k][l]] for i=1:3, j=1:3, k=1:3, l=1:3]
    
    return C4
end