module SFRPGui

using GLMakie
using SFRPMicroMechanics
S = SFRPMicroMechanics
using DataInterpolations


function start_gui(; N = 360, mandel = true)

    fig = Figure()
    ax = PolarAxis(fig[1,1])
    
    BOTTOM_BAR = GridLayout(fig[2,1], tellwidth = false)

    BOTTOM_BAR[1,1] = sg = SliderGrid(
        fig,
        (label = "Em", range = (0.1:0.1:5.0), format = "{:.1f}", startvalue = 3.0),
        (label = "nu_m", range = (0.01:0.01:0.49), format = "{:.3f}", startvalue = 0.3),
        (label = "Ef", range = (20.0:0.5:100.0), format = "{:.1f}", startvalue = 70.0),
        (label = "nu_f", range = (0.01:0.001:0.49), format = "{:.3f}", startvalue = 0.2),
        (label = "vf", range = (0.01:0.01:0.99), format = "{:.2f}", startvalue = 0.1),
        (label = "AR", range = (0.5:0.5:100), format = "{:.1f}", startvalue = 20.0),
        (label = "A11", range = (0.33:0.01:1.0), format = "{:.2f}", startvalue = 0.7),
        (label = "A22", range = (0.0:0.01:0.5), format = "{:.2f}", startvalue = 0.25),

    )

    sliderobservables = [s.value for s in sg.sliders]
    
    angles = LinRange(0, 360, N)
    Emods = Observable(zeros(size(angles)))
    on(sliderobservables[7]) do a11
         
        a22_max = min(1 - a11, a11)
        a22_min = round((1-a11) / 2; digits= 2)
        a22 = sliderobservables[8][]
        a22 = clamp(a22_min, a22_max, a22)
         
        sg.sliders[8].range = a22_min:0.01:a22_max
        sliderobservables[8][] = set_close_to!(sg.sliders[2], a22)
    end
    onany(sliderobservables...) do slvals...
        Em = slvals[1]
        nu_m = slvals[2]
        Ef = slvals[3]
        nu_f = slvals[4]
        vf = slvals[5]
        AR = slvals[6]
        a11 = slvals[7]
        a22 =slvals[8]
        Emods[] = compute_emod(angles, Em, nu_m, Ef, nu_f, vf, AR, a11, a22;mandel)
    end

    lines!(ax, deg2rad.(angles), Emods)

    fig
end

function compute_emod(phi, Em, nu_m, Ef, nu_f, vf, AR, A11, A22; 
                    
                    mandel = true,
                    )
    pm = S.IsotropicElasticParameters(Em, nu_m)
    pf = S.IsotropicElasticParameters(Ef, nu_f)
    
    # Cm = S.stiffness_matrix_voigt(pm; mandel)
    # Cf = S.stiffness_matrix_voigt(pf; mandel)

    fibers = [S.FiberPhase(pf, vf, AR, S.SpheroidalInclusion(), nothing)]

    Cmt = S.mori_tanaka(pm, fibers;mandel = true, symmetrize = true)
    
    a2= S.OrientationTensor(A11, A22)
    CT = S.HybridClosure
    Cavg = S.orientation_average(Cmt, a2; closure_type = CT)
    pavg = S.extract_orthotropic_constants(Cavg)
    return [S.apparent_modulus(pavg, ang) for ang in phi]
end


end #module