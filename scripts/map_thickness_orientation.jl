

function map_orientation_tensor(ttot, mesh; 
                    T = Float32,
                    nip = 5, 
                    pid = nothing, 
                    nplane = 1,
                    nthick = 5,
                    fname = "mapped.k", 
                    mat_215 = true,
                    dir = nothing)
    #ttot is a function f(t)-> (a11(t), a22(t))
    #mesh is a mesh file containing nodes, shell elements, section, part cards
    #produces a .k file with *initial_stress_shell for each element

    mesh_data = read_ls_dyna(mesh)
    dir = isnothing(dir) ? nothing : normalize(dir)

    #precalc 

    tgauss = T.(gauss_pts(nthick))

    hisv_dict=  Dict{Int, Tuple{T, T, T}}()
    for (i,ti) in enumerate(tgauss)
        #thickness pos, a11, a22
        hisvals = zeros(T, 18)
        (a11, a22) = ttot(ti)
        hisvals[9] = a11
        hisvals[10] = a22
        
        push!(hisv_dict, i => (ti, hisvals))
    end


    if mat_215
        hisv_dict = Dict{Int, T}(

        )
    end
    open(fname, "w+") do io

        println(io, "*INITIAL_STRESS_SHELL")
        for element in mesh_data.elements
            nodes = mesh_data.nodes[element.nodes...]
            e21 = nodes[2] - nodes[1]
            e41 = nodes[4] - nodes[1]
            normal = normalize(cross(e21, e41))
            if !isnothing(dir)
                xdir = normalize(e21)
                projdir = normalize(dir - dot(dir, normal)*dir)
                cosa = dot(xdir, projdir)
                sina = norm(cross(xdir, projdir))
                
            else
                beta = 0.0
            end
            # eid, nplane, nthick, nhisv, ntensr, large, nthint, nthhsv

            println(io, lpad.(["$eid", "$nplane", "$nhisv", "0", "0", "0", "0"], 10))
            
            print_initial_stress_shell_for_element(io, hisv_dict; eid, nplane, nthick, )

        end
    end


end


function print_initial_stress_shell_for_element(io::IO, hisv_dict; 
                    T = Float32, 
                    nplane = 1, 
                    nthick = 3,  
                    )
    
    if !isempty(hisv_dict)
        nhisv = isempty(hisv_dict) ? 0 : maximum(keys(hisv_dict))
        nhisv_lines = ceil(Int, nhisv / 8)
        for i in 1:nhisv_lines * 8
            if !haskey(hisv_dict, i)
                push!(hisv_dict, i => zero(T))
            end
        end
    end
    
    
    
    println(io, lpad.(("$eid", "$nplane", "$nthick", "$nhisv", "0", "0", "0", "0"), 10))
    for i in 1:nthick
        for _ in 1:nplane
            println(io, lpad.( ("$(tgauss[i])", ), 10)) #no sigxx
            for j in 1:nhisv_lines
                println(io, lpad.(["$(hisv_dict[k])" for k in (j-1)*8+1 : j*8], 10))
            end
        end
    end
    return nothing
end



"""
    gauss_pts(n=3)

Returns a tuple of Gauss-Legendre integration points for shell thickness.
Calculates points analytically for n <= 5 or via Newton-Raphson for n > 5.
"""
function gauss_pts(n::Int=3)
    if n == 1
        return (0.0,)
    elseif n == 2
        val = 1.0 / sqrt(3.0)
        return (-val, val)
    elseif n == 3
        val = sqrt(0.6)
        return (-val, 0.0, val)
    elseif n == 4
        v1 = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0)
        v2 = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0)
        return (-v2, -v1, v1, v2)
    elseif n == 5
        v1 = 1/3 * sqrt(5 - 2 * sqrt(10/7))
        v2 = 1/3 * sqrt(5 + 2 * sqrt(10/7))
        return (-v2, -v1, 0.0, v1, v2)
    else
        # Fallback for high n: Numerical calculation
        return Tuple(_compute_gauss_roots(n))
    end
end

function _compute_gauss_roots(n::Int)
    roots = zeros(n)
    m = (n + 1) ÷ 2
    for i in 1:m
        # Initial guess for the root
        z = cos(pi * (i - 0.25) / (n + 0.5))
        z_old = 0.0
        
        # Newton-Raphson iteration
        while abs(z - z_old) > 1e-15
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
            end
            # p1 is the Legendre polynomial, pp is its derivative
            pp = n * (z * p1 - p2) / (z^2 - 1.0)
            z_old = z
            z = z_old - p1 / pp
        end
        roots[i] = -z
        roots[n + 1 - i] = z
    end
    return roots
end


