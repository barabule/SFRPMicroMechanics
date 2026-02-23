function linear_closure(a2)
    
    SArray{Tuple{3,3,3,3}}( -1/35 * (δ(i, j) * δ(k, l) + δ(i,k) * δ(j, l) + δ(i, l) * δ(j, k) +
                            1/7 * (δ(i, j) * a2[k, l] + 
                                   δ(k, l) * a2[i, j] +
                                   δ(i, k) * a2[j, l] +
                                   δ(j, l) * a2[i, k] +
                                   δ(i, l) * a2[j, k] + 
                                   δ(j, k) * a2[i, l]))
                        for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

end

function quadratic_closure(a2)
    return SArray{Tuple{3,3,3,3}}(a2[i,j] * a2[k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end

"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
function hybrid_closure(a2)

    f = 27 * det(a2)

    a4_l = linear_closure(a2)
    a4_q = quadratic_closure(a2)

    return f * a4_l + (1 - f) * a4_q
end


function HL1_closure(a2)
    b(i, j) = sum(a2[i, m] * a2[m, j] for m in 1:3)

    SArray{Tuple{3,3,3,3}}(2/5 * (δ(i, j) * a2[k, l] + δ(k, l) * a2[i, j]) -
                           1/5 * (a2[i, j] * a2[k, l]) +
                           3/5 * (a2[i, k] * a2[j, l] + a2[i, l] * a2[j, k]) -
                           2/5 * (δ(i, j) * b(k, l) + δ(k, l) * b(i, j))
                           for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end


function HL2_closure(a2)
    
    b(i, j) = sum(a2[i, m] * a2[m, j] for m in 1:3)
    
    e(i, j, k, l) = exp(2 * (1 - 3 * a2[i, j] * a2[k, l]) / (1 - a2[i, j] * a2[k, l]))

    SArray{Tuple{3,3,3,3}}(26/315 * (δ(i,j)*δ(k,l) + δ(i,k) * δ(j, l) + δ(i, l) * δ(j, k)) *  e(i, j, k, l) +
                            16/63 * (a2[i, j] * δ(k, l) + a2[k, l] * δ(i,j)) * e(i, j, k, l) -
                            4/21  * (a2[i, k] * δ(j, l) + a2[j,l] * δ(i, k) + a2[i, l] * δ(j, k) + a2[j, k]* δ(i, l)) * 
                                     e(i, j, k, l) +
                            (a2[i,j]* a2[k,l] + a2[i,k] * a2[j, l] + a2[i, l] * a2[j, k]) - 
                            2 / (δ(i,j) * b(k,l) + δ(k,l) * b(i,j)) * b(i,j) * (b(k,l))
                                for i in 1:3, j in 1:3, k in 1:3, l in 1:3    
                                )
end