
abstract type InclusionGeometry end

struct SphericalInclusion <:InclusionGeometry
end

SphericalInclusion(nu, AR) = SphericalInclusion(nu)

struct SpheroidalInclusion <:InclusionGeometry #ellypsoidal
end

struct NeedleInclusion <: InclusionGeometry #thin cylinder with effectively infinite aspect ratio
end



struct DiscInclusion <: InclusionGeometry # thick disc / cylinder
end

struct ThinDiscInclusion <: InclusionGeometry #infinitely thin cylinder
end


function eshelby_tensor(geom::SpheroidalInclusion, nu::T1, AR::T2) where {T1<:Real, T2<:Real}
    T = promote_type(T1, T2)
    a = AR
    if a ≈ 1 
        return eshelby_tensor(SphericalInclusion(), nu, AR)
    end
    a2 = a * a
    
    if a > 1 #prolate sphere
        g = a / sqrt((a2 - 1)^3) * (a * sqrt(a2 - 1) - acosh(a))
    else #oblate 
        g = a / sqrt((1 - a2)^3) * (acos(a) - a * sqrt(1 - a2))
    end
    

    S = MArray{Tuple{3,3,3,3}, T}(zeros(3,3,3,3))

    S[1,1,1,1] = 1 / (2 * (1 - nu)) * (1 - 2nu + (3a2 - 1)/ (a2 - 1) - (1 - 2nu + 3a2 / (a2 - 1)) * g)

    S[2,2,2,2] = S[3,3,3,3] = 3 / (8 * (1 - nu)) * a2 / (a2 - 1) + 1/ (4 * (1 - nu)) * (1 - 2nu - 9/(4 * (a2 - 1))) * g

    S[2,2,3,3] = S[3,3,2,2] = 1 / (4 * (1 - nu)) * (a2 / (2 * (a2 - 1)) - (1 - 2nu + 3/(4 * (a2 -1))) * g)

    S[2,2,1,1] = S[3,3,1,1] = -1 / (2 * (1 - nu)) * a2 / (a2 - 1) + 1/(4 * (1 - nu)) * (3a2 / (a2 - 1) - (1 - 2nu)) * g
    


    S[1,1,2,2] = S[1,1,3,3] = -1 / (2 * (1 - nu))  * (1 - 2nu + 1 / (a2 - 1)) +
                               1 / (2 * (1 - nu)) * (1 -2nu + 3 / (2 * (a2 - 1))) * g 

    S[2,3,2,3] = S[3,2,3,2] = S[3,2,2,3] = S[2,3,3,2] =  1 / (4 * (1 - nu)) * (a2 / (2 * (a2 - 1)) + (1 - 2nu - 3/(4 * (a2 -1 ))) * g)
    

    S[1,2,1,2] = S[1,3,1,3] = S[2,1,2,1] = S[2,1,1,2] = S[1,2,2,1] = S[3,1,3,1] = S[3,1,1,3] = S[1,3,3,1] = 
                1 / (4 * (1 - nu)) * (1 -2nu - (a2 + 1)/ (a2 - 1) - 1/2 * (1 - 2nu  - 3 * (a2 + 1) / (a2 - 1)) * g)
    


    return SArray{Tuple{3,3,3,3}, T}(S)
end





function eshelby_tensor(geom::SphericalInclusion, nu::T, AR=nothing) where T<:Real
    
    fac = inv(15 * (1 - nu))
    
    S = MArray{Tuple{3,3,3,3}, T}(zeros(3,3,3,3))
    S[1,1,1,1] = S[2,2,2,2] = S[3,3,3,3] = (7 - 5nu) * fac
    S[1,1,2,2] = S[1,1,3,3] = S[2,2,1,1] = S[2,2,3,3] = S[3,3,1,1] = S[3,3,2,2] = (5nu - 1) * fac
    S[1,2,1,2] = S[1,2,2,1] = S[2,1,1,2] = S[2,1,2,1] =
    S[2,3,2,3] = S[2,3,3,2] = S[3,2,2,3] = S[3,2,3,2] = 
    S[1,3,1,3] = S[1,3,3,1] = S[3,1,1,3] = S[3,1,3,1] = (4 - 5nu) *fac

    return SArray{Tuple{3,3,3,3}, T}(S)
    # return SArray{Tuple{3,3,3,3}}(fac * ((5nu - 1) * (δ(i, j) * δ(k, l)) +
    #                                      (4 -5nu) * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    #                                     )
    #                             for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

end


function eshelby_tensor(geom::DiscInclusion, nu::T1, AR::T2) where {T1<:Real, T2<:Real}
    a = AR
    T = promote_type(T1, T2)
    S = MArray{Tuple{3,3,3,3}, T}(zeros(3,3,3,3))

    S[1,1,1,1] = 1 - (1 - 2nu) / (4 * (1 - nu)) * pi * a
    S[2,2,2,2] = S[3,3,3,3] = -(13 - 8nu) / (32 * (1- nu)) * pi * a
    S[2,2,3,3] = S[3,3,2,2] = (8nu - 1) / (32 * (1 - nu)) * pi * a
    S[2,2,1,1] = S[3,3,1,1] = (2nu - 1) / (8 * (1 - nu)) * pi * a
    S[1,1,2,2] = S[1,1,3,3] = nu / (1 - nu) * (1 - (1 + 4nu)/ (8nu) * pi * a)
    S[2,3,2,3] = (7 - 8nu) / (32 * (1 - nu)) * pi * a
    S[1,2,1,2] = S[1,3,1,3] = 1/2 * (1 - (2 - nu) / (4 * (1 - nu)) *  pi * a)

    return SArray{Tuple{3,3,3,3}, T}(S)
end

function eshelby_tensor(geom::ThinDiscInclusion, nu::T, AR=nothing) where T<:Real
    
    S = MArray{Tuple{3,3,3,3}, T}(zeros(3,3,3,3))
    S[1,1,1,1] = 1
    S[1,1,2,2] = S[1,1,3,3] = nu / (1 - nu)
    S[1,2,1,2] = S[1,3,1,3] = 1/2

    return SArray{Tuple{3,3,3,3}, T}(S)
end


function eshelby_tensor(geom::NeedleInclusion, nu::T, AR=nothing) where T<:Real
    
    S = MArray{Tuple{3,3,3,3}, T}(zeros(3,3,3,3))

    S[1,1,1,1] = (1 - 2nu)/ (4 * (1 - nu))
    S[2,2,2,2] = S[3,3,3,3] = (5 - 4nu) / (8 * (1 - nu))
    S[2,2,3,3] = S[3,3,2,2] = (4nu - 1)/ (8 * (1 - nu))
    S[2,2,1,1] = S[3,3,1,1] = nu / (2 * (1 - nu))
    S[2,3,2,3] = (3 - 4nu) / (8 * (1 - nu))
    S[1,2,1,2] = S[1,3,1,3] = 1/4

    return SArray{Tuple{3,3,3,3}, T}(S)
end