abstract type AbstractClosure end #all 4th order closures


function is_symmetric(A4; tol=1e-9)
    # Check all permutations of i, j, k, l
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        val = A4[i,j,k,l]
        # Test just a few key permutations to ensure full symmetry
        perms = [(j,i,k,l), (i,j,l,k), (k,l,i,j), (l,k,j,i)]
        for p in perms
            if abs(val - A4[p...]) > tol
                return false
            end
        end
    end
    return true
end

function is_normalized(A4, a2; tol=1e-9)
    # contracted A4 == a2
    # Contract A4: A_ijkk
    a2_recovered = zeros(3, 3)
    for i in 1:3, j in 1:3
        a2_recovered[i,j] = sum(A4[i,j,k,k] for k in 1:3)
    end
    
    # Check if recovered a2 matches input a2 and if trace(a2) == 1
    is_consistent = norm(a2_recovered - a2) < tol
    is_normalized = abs(tr(a2_recovered) - 1.0) < tol
    
    return is_consistent && is_normalized
end

function is_positive_semidef(A4; tol=1e-9)
    # Check diagonal components (must be non-negative)
    for i in 1:3, j in 1:3
        if A4[i,i,j,j] < -tol
            return false
        end
    end
    A66= convert_3333_to_66(A4;mandel = false)
    eigv = eigvals(A66)
    any(eigv .< -tol) && return false 
    return true
end



"""
    test_invariance(closure_func, a2; tol=1e-8)

Tests if a closure approximation is frame-invariant.
Returns true if Closure(R*a2*R') == Rotate(Closure(a2)).
"""
function is_invariant(CT::Type{<:AbstractClosure}, a2; tol=1e-8)
    # 1. Generate a random 3x3 Rotation Matrix (R)
    # We use QR decomposition of a random matrix to get an orthogonal matrix
    Q, R_mat = qr(randn(3, 3))
    R = Matrix(Q)
    if det(R) < 0  # Ensure it's a rotation, not a reflection
        R[:, 1] .*= -1
    end

    # 2. Path A: Rotate a2 first, then apply closure
    a2_rotated = R * a2 * R'
    A4_from_rotated_a2 = CT(a2_rotated)

    # 3. Path B: Apply closure to original a2, then rotate the resulting A4
    A4_original = CT(a2)
    A4_rotated_manually = zeros(3, 3, 3, 3)

    # Perform the 4th-order tensor rotation: A'_{ijkl} = R_im R_jn R_kp R_lq A_mnpq
    for i=1:3, j=1:3, k=1:3, l=1:3
        val = 0.0
        for m=1:3, n=1:3, p=1:3, q=1:3
            val += R[i,m] * R[j,n] * R[k,p] * R[l,q] * A4_original[m,n,p,q]
        end
        A4_rotated_manually[i,j,k,l] = val
    end

    # 4. Compare results
    diff = norm(A4_from_rotated_a2 - A4_rotated_manually)
    
    # println("Rotation Invariance Residual: ", diff)
    return diff < tol
end



function test_closure_approximation(a::AbstractOrientationTensor, CT::Type{<:AbstractClosure}; tol = 1e-9) 
    a2 = to_matrix(a)
    A4 = closure(a, CT)

    return is_symmetric(A4;tol) && is_normalized(A4, a2; tol) && is_positive_semidef(A4;tol) && is_invariant(CT, a2)
end


abstract type AbstractSimpleClosure <: AbstractClosure end

struct LinearClosure     <: AbstractSimpleClosure end
struct QuadraticClosure  <: AbstractSimpleClosure end
struct HybridClosure     <: AbstractSimpleClosure end
struct HL1Closure        <: AbstractSimpleClosure end
struct HL2Closure        <: AbstractSimpleClosure end

function closure(a::AbstractOrientationTensor, CT::Type{<:AbstractSimpleClosure})
    a2 = to_matrix(a)

    return CT(a2)
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

function closure(a::AbstractOrientationTensor, CT::Type{<:AbstractOrthotropicClosure})
    res = decompose_eigenvalue(a) #eigenvalue decomposition and rotation matrix 3x3
    aeig = res.tensor
    R33 = res.rotation
    R66 = convert_rot_33_to_66(R33) #rotation matrix 6x6 voigt
    a11, a22 = aeig.a11, aeig.a22
    c66 = CT(a11, a22) #compute the actual closure
    
    return convert_66_to_3333(R66 * c66 * R66'; mandel = false)#finally rotate back
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
            C[m, 10] * a2 * a2 * a2  for m in 1:3)

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


### INVARIANT BASED CLOSURES

abstract type InvariantBasedClosure <: AbstractClosure end
struct IBOF<:InvariantBasedClosure end

function get_invariants(a::AbstractOrientationTensor)
    (a1, a2, a3, a4, a5, a6) = get_all_coefficients(a)

    II = 1/2 * ((a1 + a2 + a3)^2 - (a1^2 + a2^2 + a3^2 + 2a4^2 + 2a5^2 + 2a6^2))

    III = (a1 * a2 * a3 + 2a4 * a5 * a6 - a1 * a4^2 - a2 * a5^2 - a3 * a6^2)

    return (;II, III)
end

function beta_coefficients(II, III)
    
    # β3, β4, β6

    Ct =     [0.24940908165786e2      -0.497217790110754e0            0.234146291570999e2;
             -0.435101153160329e3        0.234987975114050e2           -0.412048043372534e3;
              0.372389335663877e4       -0.391044251397838e3            0.319553200392089e4;
              0.703443657916476e4        0.153965820593506e3            0.573259594331015e4;
              0.823995187366106e6        0.152772950743819e6           -0.485212803064813e5;
             -0.133931929894245e6       -0.213755248785646e4           -0.605006113515592e5;
              0.880683515327916e6       -0.400138947092812e4           -0.477173740017567e5;
             -0.991630690741981e7       -0.185949305922308e7            0.599066486689836e7;
             -0.159392396237307e5        0.296004865275814e4           -0.110656935176569e5;
              0.800970026849796e7        0.247717810054366e7           -0.460543580680696e8;
             -0.237010458689252e7        0.101013983339062e6            0.203042960322847e7;
              0.379010599355267e8        0.732341494213578e7           -0.556606156734835e8;
             -0.337010820273821e8       -0.147919027644202e8            0.567424911007837e9;
              0.322219416256417e5       -0.104092072189767e5            0.128967058686204e5;
             -0.257258805870567e9       -0.635149929624336e8           -0.152752854956514e10;
              0.214419090344474e7       -0.247435106210237e6           -0.499321746092534e7;
             -0.449275591851490e8       -0.902980378929272e7            0.132124828143333e9;
             -0.213133920223355e08       0.7249697968073995e7          -0.162359994620983e10;
              0.157076702372204e10       0.487093452892595e9            0.792526849882218e10;
             -0.232153488525298e05       0.138088690964946e5            0.466767581292985e4;
             -0.395769398304473e10      -0.160162178614234e0           -0.128050778279459e11 ]

    C = Ct'
    bs = zeros(3)
    for i in eachindex(bs)
        bs[i] = C[i, 1] + C[i, 2] * II + C[i, 3] * II^2 + C[i, 4] * III + C[i, 5] * III^2 + 
                C[i, 6] * II * III + C[i, 7] * II^2 * III + 
                C[i,8] * II * III^2 + C[i, 9] * II^3 + C[i,10] * III^3 + C[i,11] * II^3 * III + C[i, 12] * II^2 * III^2 + 
                C[i, 13] * II * III^3 + C[i, 14] * II^4 + C[i, 15] * III^4 + C[i, 16] * II^4 * III +
                C[i, 17] * II^3 * III^2 + C[i, 18] * II^2 * III^3 +
                C[i, 19] * II * III^4 + C[i, 20] * II^5 + C[i, 21] * III^5 
    end

    (β3, β4, β6) = bs

    β1 = -3/35 + 3/25 * β3 * (1/7 + 4/7 *II + 8/3*III) - 3/25 * β4 * (1 - 8/3 * II - 14/3 * III) -
          3/25 * β6 * (1/7 - 24/21 * III - 4/7 * II + 16/3 * II * III + 8/7 * II^2)

    β2 = 6/7 - 6/35 * β3 * (1 + 4II) + 6/5 * β4 * (1/6 - II) - 6/7 * β6 * (-1/5 + 2/3 * III + 4/5 * II - 8/5 * II^2)

    β5 = -4/5 * β3 - 7/5 * β4 - 6/5 * β6 * (1 - 4/3 * II)
    
    return (β1, β2, β3, β4, β5, β6)
end

function InvariantBasedOptimalFittedClosure(a::AbstractOrientationTensor)
    # a2 = to_matrix(a)

    invariants = get_invariants(a)
    (β1, β2, β3, β4, β5, β6) = beta_coefficients(invariants.II, invariants.III)
    a2 = to_matrix(a)
    b(i,j) = sum(a2[i,m] * a2[m, j] for m in 1:3)

    a4 = SArray{Tuple{3,3,3,3}}(
            1/3 * β1 * (δ(i,j)*δ(k,l) + δ(i,k) * δ(j,l) + δ(i,l) * δ(j,k)) +
            1/6 * β2 * (δ(i,j) * a2[k,l] + δ(k,l) * a2[i,j] + δ(i,k) * a2[j,l] + δ(i,l) * a2[j,k] + δ(j,k) * a2[i,l]) +
            1/3 * β3 * (a2[i,j] * a2[k,l] + a2[i,k] * a2[j,l] + a2[i,l] * a2[j,k]) +
            1/6 * β4 * (δ(i,j) * b(k, l) + δ(k,l) * b(i, j) + δ(i, k) * b(j, l) + δ(j,l) * b(i, k) + δ(i,l) * b(j,k) + δ(j,k) * b(i,l)) +
            1/6 * β5 * (a2[i,j] * b(k,l) + a2[k,l] * b(i,j) + a2[i,k] * b(j,l) + a2[j,l] * b(i,k) + a2[i,l] * b(j,k) + a2[j,k] * b(i,l)) +
            1/6 * β6 * (b(i,j) * b(k,l) + b(i,k)*b(j,l) + b(i,l)*b(j,k))
            for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    )
    return a4
end


IBOF(a) = InvariantBasedOptimalFittedClosure(a)

function closure(a::AbstractOrientationTensor, CT::Type{<:InvariantBasedClosure})
    return CT(a)
end