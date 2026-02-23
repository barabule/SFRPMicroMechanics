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


function smooth_orthotropic_closure(a::OrientationTensor)
    a1, a2 = a.a11, a.a22

    Ct =       @SMatrix [-0.15   -0.15    0.60;
                         1.15    1.15   -0.60;
                        -0.10   -0.10   -0.60]
    C = Ct'
    A11 = C[1, 1] + C[1, 2] * a1 + C[1, 3] * a2
    A22 = C[2, 1] + C[2, 2] * a1 + C[2, 3] * a2
    A33 = C[3, 1] + C[3, 2] * a1 + C[3, 3] * a2
    
    return compute_eigenvalue_closure_matrix(a1, a2, A11, A22, A33)
end

ORS_closure(a::OrientationTensor) = smooth_orthotropic_closure(a)

function fitted_orthotropic_closure(a::OrientationTensor)
    a1, a2 = a.a11, a.a22

    Ct = @SMatrix [ 0.060964    0.124711   1.228982;
                   0.371243   -0.389402  -2.054116;
                   0.555301    0.258844   0.821548;
                  -0.369160    0.086169  -2.260574;
                   0.318266    0.796080   1.053907;
                   0.371218    0.544992   1.819756]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end

ORF_closure(a::OrientationTensor) = fitted_orthotropic_closure(a)

function low_interaction_orthotropic_closure(a::OrientationTensor)
    a1, a2 = a.a11, a.a22

    Ct = @SMatrix [ 0.104753    0.162210     1.288896;
                   0.346874   -0.451257    -2.187810;
                   0.544970    0.286639     0.899635;
                  -0.563168   -0.028702    -2.402857;
                   0.471144    0.864008     1.133547;
                   0.491202    0.652712     1.975826]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end

ORL_closure(a::OrientationTensor) = low_interaction_orthotropic_closure(a)


function wide_range_interaction_orthotropic_closure(a::OrientationTensor)

    a1, a2 = a.a11, a.a22

    Ct = @SMatrix [0.070055    0.115177        1.249811;
                   0.339376   -0.368267       -2.148297;
                   0.590331    0.252880        0.898521;
                  -0.396796    0.094820       -2.290157;
                   0.333692    0.800181        1.044147;
                   0.411944    0.535224        1.934914]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end


ORW_closure(a::OrientationTensor) = wide_range_interaction_orthotropic_closure(a)


function ortho_poly_6(a1, a2, C)
    Aii = (C[m, 1] + C[m, 2] * a1 + C[m, 3] * a1 * a1 +
                     C[m, 4] * a2 + C[m, 5] * a2 * a2 +
                                    C[m, 6] * a1 * a2 for m in 1:3) 

    return compute_eigenvalue_closure_matrix(a1, a2, Aii...)
end


function compute_eigenvalue_closure_matrix(a1, a2, A11, A22, A33)

    A23 = A44 = 1/2 * (1 - 2a1 + A11 - A22 - A33)
    A13 = A55 = 1/2 * (1 - 2a2 - A11 + A22 - A33)
    A12 = A66 = 1/2 * (-1 + 2a1 + 2a2 - A11 - A22 + A33)

    return  @SMatrix [A11 A12 A13  0   0   0;
                      A12 A22 A23  0   0   0;
                      A13 A23 A33  0   0   0;
                       0   0   0  A44  0   0;
                       0   0   0   0  A55  0;
                       0   0   0   0   0  A66]
end