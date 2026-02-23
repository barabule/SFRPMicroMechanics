function linear_closure(a)
    T = eltype(a)
    a4 = MArray{Tuple{3,3,3,3}}(zeros(T, 3, 3, 3, 3))

    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        # Linear Component
        term1 = (δ(i,j)*δ(k,l) + δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) / 24.0

        term2 = (a[i,j]*δ(k,l) + a[i,k]*δ(j,l) + a[i,l]*δ(j,k) + 
                 a[k,l]*δ(i,j) + a[j,l]*δ(i,k) + a[j,k]*δ(i,l)) / 6.0
        
        a4[i,j,k,l] = -term1 + term2
    end
    return SArray{Tuple{3,3,3,3}}(a4)
end

function quadratic_closure(a2)
    T = eltype(a2)
    a4 = MArray{Tuple{3,3,3,3}}(zeros(T, 3, 3, 3, 3))

    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        a4[i, j, k, l] = a2[i, j] * a2[k, l]
    end
    return SArray{Tuple{3,3,3,3}}(a4)
end

"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
function hybrid_closure(a2)
    f = -0.5
    for i in 1:3
        for j in 1:3
            f += 1.5 * a2[i, j] * a2[j, i]
        end
    end
    a4_l = linear_closure(a2)
    a4_q = quadratic_closure(a2)

    return (1 - f) * a4_l + f * a4_q
end