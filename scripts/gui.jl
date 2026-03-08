module SFRPGui

using GLMakie
using SFRPMicroMechanics
S = SFRPMicroMechanics
using DataInterpolations


function start_gui(; 
            N = 360, 
            mandel = true,
            BAR_WIDTH = 600,
            )

    fig = Figure()
    ax = PolarAxis(fig[1,1])
    
    
    BOTTOM_BAR = GridLayout(tellwidth = false, width = BAR_WIDTH)
    fig.layout[2,1] = BOTTOM_BAR


    closure_menu = Menu(fig, options=[("Linear", S.LinearClosure()), 
                                      ("Quadratic", S.QuadraticClosure()),
                                      ("Hybrid", S.HybridClosure())], 
                                        default = "Hybrid")
    matrix_sliders = SliderGrid(fig,
            (label = "Eₘ [GPa]", range = (0.5:0.1:5.0), startvalue = 3.0),
            (label = "νₘ", range=  (0.01:0.01:0.49), startvalue  =0.3),
            width = BAR_WIDTH/2
            )

    orientation_sliders = SliderGrid(fig,
            (label = "A₁₁", range = (0.33:0.01:1.0), startvalue = 0.7),
            (label = "A₂₂", range = (0.15:0.01:0.3), startvalue = 0.2),
            width = BAR_WIDTH/2
    )

    BOTTOM_BAR[1,1] = vgrid!(
        Label(fig, "Matrix", halign = :left),
        matrix_sliders,
        Label(fig, "Orientation", halign = :left),
        orientation_sliders,
        hgrid!(Label(fig, "Closure"), closure_menu),
        width = BAR_WIDTH/2,
        tellwidth = false,
    )

    

    fiber_isotropic_sliders = SliderGrid(fig,
            (label = "E", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 70.0),
            (label = "nu", range = (0.01:0.01:0.49), format = "{:.2f}", startvalue = 0.22),
            (label = "fract", range = (0.01:0.01:1.0), format = "{:.2f}", startvalue = 0.1),
            (label = "aspect", range = (0.1:0.1:50.0), format = "{:.2f}", startvalue = 10.0),
            width = BAR_WIDTH/2,
            
            )

    fiber_trans_sliders = SliderGrid(fig,
            (label = "E₁ [GPa]", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 230.0),
            (label = "E₂ [GPa]", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 35.0),
            (label = "ν₂₁", range = (0.01:0.01:0.2), format = "{:.2f}", startvalue = 0.03),
            (label = "ν₃₂", range = (0.2:0.01:0.49), format = "{:.2f}", startvalue = 0.39),
            (label = "G₁₂ [GPa]", range = (1.0:0.5:200.0), format = "{:.1f}", startvalue = 50.0),
            (label = "vol_f", range = (0.01:0.01:1.0), format = "{:.2f}", startvalue = 0.1),
            (label = "aspect", range = (0.1:0.1:50.0), format = "{:.2f}", startvalue = 10.0),
            width = BAR_WIDTH/2,
            )


    inclusion_menu = Menu(fig, options = [("Spheroidal", S.SpheroidalInclusion()),
                                          ("Sphere", S.SphericalInclusion()),
                                          ("Needle", S.NeedleInclusion()),
                                          ("Disk", S.DiscInclusion()),
                                          ("Thin Disk", S.ThinDiscInclusion())],
                                        default = "Spheroidal")

    offscreen = GridLayout(bbox = (-200, -100,0, 100))
    offscreen[1,1] = fiber_isotropic_sliders
    
    
    cb_isotropic = Checkbox(fig, checked = false)

    # fiber_elastic_sliders = Observable(fiber_isotropic_sliders)

     BOTTOM_BAR[1,2] = BOTTOM_RIGHT_BAR = GridLayout()
     
     BOTTOM_RIGHT_BAR[1,1] = vgrid!(
        Label(fig, text = "Fiber Properties", halign = :left),
        hgrid!(cb_isotropic, Label(fig, "Isotropic")),
        width = BAR_WIDTH/2,
        )

    BOTTOM_RIGHT_BAR[2,1] = fiber_trans_sliders

    BOTTOM_RIGHT_BAR[3,1] = hgrid!(Label(fig, "Inclusion geometry"), inclusion_menu)
        


    on(cb_isotropic.checked) do val
        
        BOTTOM_RIGHT_BAR[2,1] = val ? fiber_isotropic_sliders : fiber_trans_sliders
        offscreen[1,1] = val ? fiber_trans_sliders : fiber_isotropic_sliders
    end


    # sliderobservables = [s.value for s in sg.sliders]
    
    # angles = LinRange(0, 360, N)
    # phi = deg2rad.(angles)
    # Emods = Observable(zeros(size(angles)))
    # on(sliderobservables[7]) do a11
         
    #     a22_max = min(1 - a11, a11)
    #     a22_min = round((1-a11) / 2 + 0.01; digits= 2)
    #     # a22_min = round(max(1-a11, a11)/2; digits= 2)
    #     a22 = sliderobservables[8][]
    #     a22 = clamp(a22_min, a22_max, a22)
         
    #     sg.sliders[8].range = a22_min:0.01:a22_max
    #     sliderobservables[8][] = set_close_to!(sg.sliders[2], a22)
    # end


    # onany(sliderobservables...) do slvals...
    #     Em = slvals[1]
    #     nu_m = slvals[2]
    #     Ef = slvals[3]
    #     nu_f = slvals[4]
    #     vf = slvals[5]
    #     AR = slvals[6]
    #     a11 = slvals[7]
    #     a22 =slvals[8]
    #     Emods[] = compute_emod(angles, Em, nu_m, Ef, nu_f, vf, AR, a11, a22;mandel)
    # end

    # lines!(ax, phi, Emods, color = :black)

    fig
end


function slider_line(parent, label, range, startvalue)
    sl = Slider(parent, range, startvalue)
    GL = GridLayout()
    GL[1,1] = hgrid!(Label(parent, "$label", halign=:left), 
                    sl, 
                    Label(parent, lift(v-> "$(round(v, digits=2))", sl.value)))
    colsize!(GL, 2, Relative(0.6))
    return (;slider = sl, gridlayout = GL)
end


function compute_emod(phi, Em, nu_m, Ef, nu_f, vf, AR, A11, A22; 
                    mandel = true,
                    )
    pm = S.IsotropicElasticParameters(Em, nu_m)
    pf = S.IsotropicElasticParameters(Ef, nu_f)
    
    # Cm = S.stiffness_matrix_voigt(pm; mandel)
    # Cf = S.stiffness_matrix_voigt(pf; mandel)

    fibers = [S.FiberPhase(pf, vf, AR, S.SpheroidalInclusion(), nothing)]

    Cmt = S.mori_tanaka(pm, fibers;mandel, symmetrize = true)
    
    a2= S.OrientationTensor(A11, A22)
    CT = S.HybridClosure
    Cavg = S.orientation_average(Cmt, a2; closure_type = CT)

    pavg = S.extract_orthotropic_constants(Cavg)
    return [S.apparent_modulus(pavg, ang) for ang in phi]
end


end #module