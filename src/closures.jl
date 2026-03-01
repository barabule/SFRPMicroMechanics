abstract type AbstractClosure end #all 4th order closures
abstract type AbstractSimpleClosure <: AbstractClosure end

struct LinearClosure     <: AbstractSimpleClosure end
struct QuadraticClosure  <: AbstractSimpleClosure end
struct HybridClosure     <: AbstractSimpleClosure end
struct HL1Closure        <: AbstractSimpleClosure end
struct HL2Closure        <: AbstractSimpleClosure end

function compute_closure(a::AbstractOrientationTensor, closure::AbstractSimpleClosure)
    a2 = to_matrix(a)

    return closure(a2)
end




function linear_closure(a2::AbstractMatrix)
    
    SArray{Tuple{3,3,3,3}}( -1/35 * (δ(i, j) * δ(k, l) + δ(i,k) * δ(j, l) + δ(i, l) * δ(j, k)) +
                            1/7 * (δ(i, j) * a2[k, l] + 
                                   δ(k, l) * a2[i, j] +
                                   δ(i, k) * a2[j, l] +
                                   δ(j, l) * a2[i, k] +
                                   δ(i, l) * a2[j, k] + 
                                   δ(j, k) * a2[i, l])
                        for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

end

LinearClosure(a2) = linear_closure(a2)



function quadratic_closure((a2::AbstractMatrix))
    
    return SArray{Tuple{3,3,3,3}}(a2[i,j] * a2[k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end


QuadraticClosure(a2) = quadratic_closure(a2)


"""
Computes the 4th-order orientation tensor using the Hybrid Closure approximation.
Balances Linear and Quadratic closures based on the degree of alignment.
"""
function hybrid_closure((a2::AbstractMatrix))
   
    f = 27 * det(a2)

    a4_l = linear_closure(a2)
    a4_q = quadratic_closure(a2)

    return f * a4_l + (1 - f) * a4_q
end

HybridClosure(a2) = hybrid_closure(a2)


function HL1_closure((a2::AbstractMatrix))
    
    b(i, j) = sum(a2[i, m] * a2[m, j] for m in 1:3)

    SArray{Tuple{3,3,3,3}}(2/5 * (δ(i, j) * a2[k, l] + δ(k, l) * a2[i, j]) -
                           1/5 * (a2[i, j] * a2[k, l]) +
                           3/5 * (a2[i, k] * a2[j, l] + a2[i, l] * a2[j, k]) -
                           2/5 * (δ(i, j) * b(k, l) + δ(k, l) * b(i, j))
                           for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
end

HL1Closure(a2) = HL1_closure(a2)

function HL2_closure((a2::AbstractMatrix))
    
    
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

HL2Closure(a2)  = HL2_closure(a2)

#eigenvalue decomp closures
abstract type AbstractOrthotropicClosure<:AbstractClosure end
struct ORS  <: AbstractOrthotropicClosure end
struct ORF  <: AbstractOrthotropicClosure end
struct ORL  <: AbstractOrthotropicClosure end
struct ORW  <: AbstractOrthotropicClosure end
struct ORW3 <: AbstractOrthotropicClosure end
struct ORFM <: AbstractOrthotropicClosure end

function compute_closure(a::AbstractOrientationTensor, closure::AbstractOrthotropicClosure)
    (aeig, R33) = decompose_eigenvalue #eigenvalue decomposition and rotation matrix 3x3
    R66 = convert_rot_33_to_66(R33) #rotation matrix 6x6 voigt
    a11, a22 = aeig.a11, aeig.a22
    c66 = closure(a11, a22) #compute the actual closure
    
    return R66 * c66 * R66'#finally rotate back
end




function smooth_orthotropic_closure(a1::T, a2::T) where T<:Real
    

    Ct =       @SMatrix [-0.15   -0.15    0.60;
                         1.15    1.15   -0.60;
                        -0.10   -0.10   -0.60]
    C = Ct'
    A11 = C[1, 1] + C[1, 2] * a1 + C[1, 3] * a2
    A22 = C[2, 1] + C[2, 2] * a1 + C[2, 3] * a2
    A33 = C[3, 1] + C[3, 2] * a1 + C[3, 3] * a2
    
    return compute_eigenvalue_closure_matrix(a1, a2, A11, A22, A33)
end

ORS(a1, a2) = smooth_orthotropic_closure(a1, a2)

function fitted_orthotropic_closure(a1::T, a2::T) where T<:Real
    

    Ct = @SMatrix [ 0.060964    0.124711   1.228982;
                   0.371243   -0.389402  -2.054116;
                   0.555301    0.258844   0.821548;
                  -0.369160    0.086169  -2.260574;
                   0.318266    0.796080   1.053907;
                   0.371218    0.544992   1.819756]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end

ORF(a1, a2) = fitted_orthotropic_closure(a1, a2)



function low_interaction_orthotropic_closure(a1::T, a2::T) where T<:Real
    
    Ct = @SMatrix [ 0.104753    0.162210     1.288896;
                   0.346874   -0.451257    -2.187810;
                   0.544970    0.286639     0.899635;
                  -0.563168   -0.028702    -2.402857;
                   0.471144    0.864008     1.133547;
                   0.491202    0.652712     1.975826]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end

ORL(a1, a2) = low_interaction_orthotropic_closure(a1, a2)


function wide_range_interaction_orthotropic_closure(a1::T, a2::T) where T<:Real
    
    Ct = @SMatrix [0.070055    0.115177        1.249811;
                   0.339376   -0.368267       -2.148297;
                   0.590331    0.252880        0.898521;
                  -0.396796    0.094820       -2.290157;
                   0.333692    0.800181        1.044147;
                   0.411944    0.535224        1.934914]
    C = Ct'

    return ortho_poly_6(a1, a2, C)
end


ORW(a1, a2) = wide_range_interaction_orthotropic_closure(a1, a2)


function wide_range_interaction_orthotropic_closure_3rd_order(a1::T, a2::T) where T<:Real
    
    Ct = @SMatrix [ -0.1480648093    -0.2106349673     0.4868019601;
                     0.8084618453     0.9092350296     0.5776328438;
                     0.3722003446    -1.2840654776    -2.2462007509;
                     0.7765597096     1.1104441966     0.4605743789;
                    -1.3431772379     0.1260059291    -1.9088154281;
                    -1.7366749542    -2.5375632310    -4.8900459209;
                     0.8895946393     1.9988098293     4.0544348937;
                     1.7367571741     1.4863151577     3.8542602127;
                    -0.0324756095     0.5856304774     1.1817992322;
                     0.6631716575    -0.0756740034     0.9512305286]
    C = Ct'
    ortho_poly_10_terms(a1, a2, C)
end

ORW3(a1, a2) = wide_range_interaction_orthotropic_closure_3rd_order(a1, a2)




function ortho_poly_6(a1, a2, C)
    Aii = (C[m, 1] + C[m, 2] * a1 + C[m, 3] * a1 * a1 +
                     C[m, 4] * a2 + C[m, 5] * a2 * a2 +
                                    C[m, 6] * a1 * a2 for m in 1:3) 

    return compute_eigenvalue_closure_matrix(a1, a2, Aii...)
end

function ortho_poly_10_terms(a1, a2, C)
    @assert size(C) == (3, 10)

    Aii = (C[m, 1] + C[m, 2] * a1 + C[m, 3] * a1 * a1 +
                     C[m, 4] * a2 + C[m, 5] * a2 * a2 +
                                    C[m, 6] * a1 * a2 +
            C[m,  7] * a1 * a1 * a2 +
            C[m,  8] * a1 * a2 * a2 +
            C[m,  9] * a1 * a1 * a1 +
            C[m, 10] * a2 * a2 * a2  )

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


function modified_fitted_orthotropic_closure(a1, a2)
    
    return compute_modified_eigenvalue_closure_matrix(a1, a2)
    
end

ORFM(a1, a2) = modified_fitted_orthotropic_closure(a1, a2)

function compute_modified_eigenvalue_closure_matrix(a1, a2)
    A23 = A44 = 0.2 - 0.2 * a1 - 0.2 * a2
    A13 = A55 =  -0.082617 + 0.646735 * a1 - 0.559003 * a1 * a1 +
                             0.357952 * a2 - 0.288046 * a2 * a2 -
                                             0.822991 * a1 * a2
    

    A22 = 0.124711  - 0.389402 * a1 + 0.258844 * a1 * a1 +
                      0.086169 * a2 + 0.796080 * a2 * a2 +
                                    + 0.544992 * a1 * a2
    a3 = 1 - a1 - a2
    A33 = a3 - A13 - A23
    A12 =  A66 = a2 - A22 - A23
    A11 = a1 - A12 - A13

    return  @SMatrix [A11 A12 A13  0   0   0;
                      A12 A22 A23  0   0   0;
                      A13 A23 A33  0   0   0;
                       0   0   0  A44  0   0;
                       0   0   0   0  A55  0;
                       0   0   0   0   0  A66]

end


function convert_rot_33_to_66(R33)
    Q = R33 #already reordered
    Q11, Q21, Q31 = Q[1,1], Q[2,1], Q[3,1]
    Q12, Q22, Q32 = Q[1,2], Q[2,2], Q[3,2]
    Q13, Q23, Q33 = Q[1,3], Q[2,3], Q[3,3]
    
    R66 = @SMatrix [Q11 * Q11  Q12 * Q12  Q13 * Q13    2Q12 * Q13     2Q11 * Q13         2Q11 * Q12;
                    Q21 * Q21  Q22 * Q22  Q23 * Q23    2Q22 * Q23     2Q21 * Q23         2Q21 * Q22;
                    Q31 * Q31  Q32 * Q32  Q33 * Q33    2Q31 * Q33     2Q31 * Q33         2Q31 * Q32;
                    Q21 * Q31  Q22 * Q32  Q23 * Q33    Q22 * Q33 + Q23 * Q32  Q23 * Q31 + Q21 * Q33  Q21 * Q32 + Q22 * Q31;
                    Q31 * Q11  Q32 * Q12  Q33 * Q13    Q32 * Q13 + Q33 * Q12  Q33 * Q11 + Q31 * Q13  Q31 * Q12 + Q32 * Q11;
                    Q11 * Q21  Q12 * Q22  Q13 * Q23    Q12 * Q23 + Q13 * Q22  Q13 * Q21 + Q11 * Q23  Q11 * Q22 + Q12 * Q21]
    return R66
end