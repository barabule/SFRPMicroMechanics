

function qij_x(angle::T) where T
    s,c = sincosd(angle)

    return SMatrix{3,3,T}([1  0  0;
                    0  c  s;
                    0 -s  c])
end

function qij_y(angle::T) where T
    s, c = sincosd(angle)

    return SMatrix{3, 3, T}([c  0  -s;
                       0  1   0;
                       s  0   c])
end

function qij_z(angle::T) where T
    s, c = sincosd(angle)

    return SMatrix{3, 3, T}([ c  s   0;
                              -s  c   0;
                               0  0   1])
end
"""
    rotation matrix for engineering notation given some rotation tensor, qij
    note: uses convention for sigma = [s11, s22, s33, s23, s13, s12]
    input: q (3x3 rotation tensor)
    output: R_s (6x6 rotation matrix)
    """
function R_sigma(q)
    # R_s = np.array([[q[0,0]**2,q[0,1]**2,q[0,2]**2,2.*q[0,1]*q[0,2],2.*q[0,0]*q[0,2], 2.*q[0,0]*q[0,1]],
    #                    [q[1,0]**2,q[1,1]**2,q[1,2]**2,2.*q[1,1]*q[1,2],2.*q[1,0]*q[1,2], 2.*q[1,0]*q[1,1]],
    #                    [q[2,0]**2,q[2,1]**2,q[2,2]**2,2.*q[2,1]*q[2,2],2.*q[2,0]*q[2,2], 2.*q[2,0]*q[2,1]],
    #                    [q[1,0]*q[2,0], q[1,1]*q[2,1], q[1,2]*q[2,2],q[1,2]*q[2,1]+q[1,1]*q[2,2], q[1,2]*q[2,0]+q[1,0]*q[2,2], q[1,1]*q[2,0]+q[1,0]*q[2,1]],
    #                    [q[0,0]*q[2,0], q[0,1]*q[2,1], q[0,2]*q[2,2],q[0,2]*q[2,1]+q[0,1]*q[2,2], q[0,2]*q[2,0]+q[0,0]*q[2,2], q[0,1]*q[2,0]+q[0,0]*q[2,1]],
    #                    [q[0,0]*q[1,0], q[0,1]*q[1,1], q[0,2]*q[1,2],q[0,2]*q[1,1]+q[0,1]*q[1,2], q[0,2]*q[1,0]+q[0,0]*q[1,2], q[0,1]*q[1,0]+q[0,0]*q[1,1]]])
    # return R_s

    R11 = q[1,1]^2
    R12 = q[1,2]^2
    R13 = q[1,3]^2
    R14 = 2q[1,2] * q[1,3]
    R15 = 2q[1,1] * q[1,3]
    R16 = 2q[1,1] * q[1,2]

    R21 = q[2,1]^2
    R22 = q[2,2]^2
    R23 = q[2,3]^2
    R24 = 2q[2,2] * q[2,3]
    R25 = 2q[2,1] * q[2,3]
    R26 = 2q[2,1] * q[2,2]

    R31 = q[3,1]^2
    R32 = q[3,2]^2
    R33 = q[3,3]^2
    R34 = 2q[3,2] * q[3,3]
    R35 = 2q[3,1] * q[3,3]
    R36 = 2q[3,1] * q[3,2]

    R41 = q[2,1] * q[3,1]
    R42 = q[2,2] * q[3,2]
    R43 = q[2,3] * q[3,3]
    R44 = q[2,3] * q[3,2] + q[2,2] * q[3,3]
    R45 = q[2,3] * q[3,1] + q[2,1] * q[3,3]
    R46 = q[2,2] * q[3,1] + q[2,1] * q[3,1]

    R51 = q[1,1] * q[3,1]
    R52 = q[1,2] * q[3,2]
    R53 = q[1,3] * q[3,3]
    R54 = q[1,3] * q[3,2] + q[1,2] * q[3,3]
    R55 = q[1,3] * q[3,1] + q[1,1] * q[3,3]
    R56 = q[1,2] * q[3,1] + q[1,1] * q[3,2]

    R61 = q[1,1] * q[2,1]
    R62 = q[1,2] * q[2,2]
    R63 = q[1,3] * q[2,3]
    R64 = q[1,3] * q[2,2] + q[1,2] * q[2,3]
    R65 = q[1,3] * q[2,1] + q[1,1] * q[2,3]
    R66 = q[1,2] * q[2,1] + q[1,1] * q[2,2]

    return SMatrix{6,6}([R11 R12 R13 R14 R15 R16;
                         R21 R22 R23 R24 R25 R26;
                         R31 R32 R33 R34 R35 R36;
                         R41 R42 R43 R44 R45 R46;
                         R51 R52 R53 R54 R55 R56;
                         R61 R62 R63 R64 R65 R66])


end

"""
    tens2eng(tens)

    returns 4th-order tensor in engineering notation
"""
function tens2eng(tens)
    T = eltype(tens)
    A = MMatrix{6,6, T}(zeros(6, 6))

    code = ((1,1), (2,2), (3,3), (2,3), (1,3), (1, 2))
    for i in 1:6
        for j in 1:6
            A[i, j] = tens[code[i]..., code[j]...]
        end
    end
    return SMatrix{6,6}(A)
end

# function eshelby_tensor(nu_m, ar)


# end