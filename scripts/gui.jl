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

#     struct ORS  <: AbstractOrthotropicClosure end
# struct ORF  <: AbstractOrthotropicClosure end
# struct ORL  <: AbstractOrthotropicClosure end
# struct ORW  <: AbstractOrthotropicClosure end
# struct ORW3 <: AbstractOrthotropicClosure end
# struct ORFM <: AbstractOrthotropicClosure end


    closure_menu = Menu(fig, options=[("Linear", S.LinearClosure), 
                                      ("Quadratic", S.QuadraticClosure),
                                      ("Hybrid", S.HybridClosure),
                                      ("IBOF", S.IBOF),
                                      ("ORS", S.ORS),
                                      ("ORF", S.ORF),
                                      ("ORW", S.ORW),
                                      ("ORW3", S.ORW3),
                                      ("ORFM", S.ORFM),
                                      ("ORL", S.ORL),
                                      ], 
                                        default = "IBOF")
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

    USE_MT_Toggle = Toggle(fig, active = true)
    averaging_method = lift(USE_MT_Toggle.active) do val 
        val ? "Mori Tanaka" : "Halpin Tsai"
    end

    BOTTOM_BAR[1,1] = vgrid!(
        hgrid!(USE_MT_Toggle, Label(fig, averaging_method)),
        Label(fig, "Matrix", halign = :left),
        matrix_sliders,
        Label(fig, "Orientation", halign = :left),
        orientation_sliders,
        hgrid!(Label(fig, "Closure"), closure_menu),
        width = BAR_WIDTH/2,
        tellwidth = false,
    )

    matrix_observables = [s.value for s in matrix_sliders.sliders]
    matrix_properties = @lift S.IsotropicElasticParameters($(matrix_observables[1]), $(matrix_observables[2]))
    # @info "matrix", matrix_properties[]
    orientation_tensor_observables = [s.value for s in orientation_sliders.sliders]
    orientation_tensor = @lift S.OrientationTensor($(orientation_tensor_observables[1]), 
                                                    $(orientation_tensor_observables[2]))

    fiber_isotropic_sliders = SliderGrid(fig,
            (label = "E [GPa]", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 70.0),
            (label = "nu", range = (0.01:0.01:0.49), format = "{:.2f}", startvalue = 0.22),
            (label = "fract", range = (0.01:0.01:1.0), format = "{:.2f}", startvalue = 0.1),
            (label = "aspect", range = (0.1:0.1:100.0), format = "{:.2f}", startvalue = 10.0),
            width = BAR_WIDTH/2,
            
            )

    fiber_iso_observables = [s.value for s in fiber_isotropic_sliders.sliders]


    fiber_trans_sliders = SliderGrid(fig,
            (label = "E₁ [GPa]", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 230.0),
            (label = "E₂ [GPa]", range = (1.0:0.5:500.0), format = "{:.1f}", startvalue = 35.0),
            (label = "ν₂₁", range = (0.01:0.01:0.2), format = "{:.2f}", startvalue = 0.03),
            (label = "ν₃₂", range = (0.2:0.01:0.49), format = "{:.2f}", startvalue = 0.39),
            (label = "G₁₂ [GPa]", range = (1.0:0.5:200.0), format = "{:.1f}", startvalue = 50.0),
            (label = "vol_f", range = (0.01:0.01:1.0), format = "{:.2f}", startvalue = 0.1),
            (label = "aspect", range = (0.1:0.1:100.0), format = "{:.2f}", startvalue = 10.0),
            width = BAR_WIDTH/2,
            )

    fiber_trans_observables = [s.value for s in fiber_trans_sliders.sliders]

    inclusion_menu = Menu(fig, options = [("Spheroidal", S.SpheroidalInclusion()),
                                          ("Sphere", S.SphericalInclusion()),
                                          ("Needle", S.NeedleInclusion()),
                                          ("Disk", S.DiscInclusion()),
                                          ("Thin Disk", S.ThinDiscInclusion())],
                                        default = "Spheroidal")

    offscreen = GridLayout(bbox = (-200, -100,0, 100))
    offscreen[1,1] = fiber_isotropic_sliders
    
    offscreen2 = GridLayout(bbox = (-200, -100,0, 100))
    
    
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
        # isotropic = val || USE_MT_Toggle.active[] 
        # if USE_MT_Toggle.active[]# keep iso always when H-T

        # end
        BOTTOM_RIGHT_BAR[2,1] = val ? fiber_isotropic_sliders : fiber_trans_sliders
        offscreen[1,1] = val ? fiber_trans_sliders : fiber_isotropic_sliders
    end


   
    
    vol_frac = Observable(0.1)
    aspect_ratio = Observable(10.0)



    angles = LinRange(0, 360, N)
    phi = deg2rad.(angles)
    

    fiber_elastic_props = Observable{S.AbstractElasticParameters}(S.IsotropicElasticParameters(70.0, 0.22))

    onany(fiber_iso_observables...) do vals...
        Ef = vals[1]
        nuf = vals[2]
        vf = vals[3]
        ar = vals[4]
        vol_frac[] = vf
        aspect_ratio[] = ar
        fiber_elastic_props[] = S.IsotropicElasticParameters(Ef, nuf)
    end

    onany(fiber_trans_observables...) do vals...
        E1 = vals[1]
        E2 = vals[2]
        nu21 = vals[3]
        nu23 = vals[4]
        G12 = vals[5]
        vf = vals[6]
        ar = vals[7]
        G23 = E2/ (2*(1+nu23))
        G31 = G12
        nu31 = nu21
        fiber_elastic_props[] = S.OrthotropicElasticParameters(;E1, E2, E3=E2,nu21, nu23, nu31, G12, G23, G31) 
        vol_frac[] = vf
        aspect_ratio[] = ar
    end

    on(orientation_sliders.sliders[1].value) do a11
         
        a22_max = min(1 - a11, a11)
        a22_min = round((1-a11) / 2 + 0.01; digits= 2)
        # a22_min = round(max(1-a11, a11)/2; digits= 2)
        # a22 = orientation_tensor_observables[2][]
        # a22 = clamp(a22_min, a22_max, a22)
        a22 = 1/2 * (a22_min + a22_max)
        a22_slider = orientation_sliders.sliders[2]
        a22_slider.range = a22_min:0.01:a22_max
        orientation_tensor_observables[2][] = set_close_to!(a22_slider, a22)
    end


    
    Emods = @lift compute_emod($matrix_properties, 
                               $fiber_elastic_props, 
                               $vol_frac,
                               $aspect_ratio,
                               $orientation_tensor,
                               $(inclusion_menu.selection),
                               $(closure_menu.selection),
                               angles;
                               mandel = true,
                               symmetrize = true,
                               mori_tanaka = $(USE_MT_Toggle.active),
                               )
    lines!(ax, phi, Emods, color = :black)


    set_close_to!(matrix_sliders.sliders[1], 3.0)

    fig
end


function compute_emod(pm::S.IsotropicElasticParameters, 
                      pf::S.AbstractElasticParameters, 
                      volume_fraction, 
                      aspect_ratio, 
                      a2::S.OrientationTensor,
                      inclusion::S.InclusionGeometry,
                      closure::Type{<:S.AbstractClosure},
                      angles::AbstractVector{<:Real};
                      mandel = true,
                      symmetrize = true,
                      mori_tanaka = true)

    # @info "pm", pm
    # @info "pf", pf
    # @info "vf", volume_fraction
    # @info "ar", aspect_ratio
    # @info "a2", a2
    # @info "inclusion", inclusion
    # @info "closure", closure

    
    if mori_tanaka
        fibers = [S.FiberPhase(pf, volume_fraction, aspect_ratio, inclusion)]
        Cmt = S.mori_tanaka(pm, fibers;mandel, symmetrize)
    else
        
        Cmt = S.halpin_tsai(pm, pf, volume_fraction, aspect_ratio)
    end
    
    Cavg = S.orientation_average(Cmt, a2; closure_type = closure, mandel)
    # @info "Cavg" 
    # display(Cavg)

    pavg = S.extract_orthotropic_constants(Cavg; mandel)
    # @info "Extracted props", pavg
    return [S.apparent_modulus(pavg, ang) for ang in angles]

end

end #module