module SFRPGui

using GLMakie
using SFRPMicroMechanics
S = SFRPMicroMechanics
using DataInterpolations


function start_gui(; 
            N = 360, 
            mandel = true,
            BAR_WIDTH = 600,
            ref_data = nothing,
            )

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Direction [°]", ylabel = "E Modulus [GPa]",
                limits= (nothing, nothing, 0, nothing))
    
    
    BOTTOM_BAR = GridLayout(tellwidth = false, width = BAR_WIDTH)
    fig.layout[2,1] = BOTTOM_BAR


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
            (label = "αₘ", range = (10:1:200), startvalue = 80),
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
            (label = "αf", range = (-5:1:50), startvalue = 5),
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
            (label = "α₁", range = (-5:1:50), startvalue = 5),
            (label = "α₂", range = (-5:1:50), startvalue = 5),
            (label = "α₃", range = (-5:1:50), startvalue = 5),
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
        
    cte_matrix = @lift S.ThermalExpansion($(matrix_observables[3])* 1e-6)
    cte_fiber  = Observable(S.ThermalExpansion( 5e-6))

    on(cb_isotropic.checked) do val
        # isotropic = val || USE_MT_Toggle.active[] 
        # if USE_MT_Toggle.active[]# keep iso always when H-T

        # end
        BOTTOM_RIGHT_BAR[2,1] = val ? fiber_isotropic_sliders : fiber_trans_sliders
        offscreen[1,1] = val ? fiber_trans_sliders : fiber_isotropic_sliders
    end


   
    
    vol_frac = Observable(0.1)
    aspect_ratio = Observable(10.0)



    angles = LinRange(0, 90, N)
    phi = deg2rad.(angles)
    

    fiber_elastic_props = Observable{S.AbstractElasticParameters}(S.IsotropicElasticParameters(70.0, 0.22))

    onany(fiber_iso_observables...) do vals...
        Ef = vals[1]
        nuf = vals[2]
        alfa_m = vals[3]
        vf = vals[4]
        ar = vals[5]
        vol_frac[] = vf
        aspect_ratio[] = ar
        fiber_elastic_props[] = S.IsotropicElasticParameters(Ef, nuf)
        cte_fiber[] = S.ThermalExpansion(alfa_m * 1e-6)
    end

    onany(fiber_trans_observables...) do vals...
        E1 = vals[1]
        E2 = vals[2]
        nu21 = vals[3]
        nu23 = vals[4]
        G12 = vals[5]
        alfa1 = vals[6]
        alfa2 = vals[7]
        alfa3 = vals[8]
        vf = vals[9]
        ar = vals[10]
        G23 = E2/ (2*(1+nu23))
        G31 = G12
        nu31 = nu21
        fiber_elastic_props[] = S.TransverseIsotropicElasticParameters(;E1, E2, nu21, nu23, nu31, G12) 
        vol_frac[] = vf
        aspect_ratio[] = ar
        cte_fiber[] = S.ThermalExpansion(alfa1 * 1e-6, alfa2 * 1e-6, alfa3 * 1e-6)
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

    
    
    res = @lift compute_emod($matrix_properties, 
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

    Emod = @lift $res.moduli
    lines!(ax, angles, Emod, color = :black)

    if !isnothing(ref_data)
        angs = ref_data.angles
        vmin = ref_data.min
        vmax = ref_data.max
        mids = @. 1/2 * (vmin + vmax)
        err = @. 1/2  * (vmax - vmin)

        errorbars!(ax, angs, mids, err, whiskerwidth = 10, color = :red, linewidth = 4)
    end


    Em = lift(matrix_properties) do pm
        pm.E_modulus
    end
    # hlines!(ax, [Em], color = :black, linestyle = :dash)
    
    cte_eff = @lift compute_effective_thermal_expansion($matrix_properties,
                                                        $fiber_elastic_props,
                                                        $cte_matrix,
                                                        $cte_fiber,
                                                        $vol_frac,
                                                        $aspect_ratio,
                                                        $orientation_tensor,
                                                        $(inclusion_menu.selection);
                                                        mandel
                                                        )




    set_close_to!(matrix_sliders.sliders[1], 3.0)

    # ptext = @lift display($(res.properties))
    # stats_text = @lift display($res.properties) * cte_to_text($cte_eff)
    stats_text = @lift sprint(show,"text/plain", $res.properties) * "\n" * cte_to_text($cte_eff)


    text!(
        ax,
        Point2f(10, 50), # 2D position in pixels
        text = stats_text,
        fontsize = 16,
        color = :black,
        space = :pixel, # Renders text in screen space
        # align = (:right, :top),
        # clip = false,
    )


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
    return (;moduli = [S.apparent_modulus(pavg, ang) for ang in angles], properties = pavg)

end

function compute_effective_thermal_expansion(pm::S.IsotropicElasticParameters, 
                                             pf::S.AbstractElasticParameters, 
                                             cte_m::S.ThermalExpansion,
                                             cte_f::S.ThermalExpansion,
                                             volume_fraction,
                                             aspect_ratio,
                                             a2::S.OrientationTensor,
                                             inclusion::S.InclusionGeometry;
                                             mandel = true,
                                             average = true,
                                             )

    fibers = [S.FiberPhase(pf, volume_fraction, aspect_ratio, inclusion)]
    ctes = [cte_m, cte_f]

    all_ctes  = S.compute_all_thermal_expansions(pm, fibers, ctes, a2; mandel, average) 
    # cte_eff = S.ThermalExpansion(pm, pf, cte_m, cte_f,volume_fraction, aspect_ratio, a2, inclusion)
    return all_ctes
end


function cte_to_text(cte_eff)

    alpha_turner = round(cte_eff.turner.alpha1 * 1e6, digits = 1)
    
    alpha_rom = round(cte_eff.rom.alpha1 * 1e6, digits = 1)
    
    alpha_kerner = round(cte_eff.kerner.alpha1 * 1e6, digits = 1)

    alpha1_shapery = round(cte_eff.shapery.alpha1 * 1e6, digits = 1)
    alpha2_shapery = round(cte_eff.shapery.alpha2 * 1e6, digits = 1)

    alpha1_chow = round(cte_eff.chow.alpha1 * 1e6, digits = 1)
    alpha2_chow = round(cte_eff.chow.alpha2 * 1e6, digits = 1)
    alpha3_chow = round(cte_eff.chow.alpha3 * 1e6, digits = 1)

    alpha1_mt = round(cte_eff.moritanaka.alpha1 * 1e6, digits = 1)
    alpha2_mt = round(cte_eff.moritanaka.alpha2 * 1e6, digits = 1)
    alpha3_mt = round(cte_eff.moritanaka.alpha3 * 1e6, digits = 1)

    io = IOBuffer()
    println(io, "Effective CTEs:")
    println(io, "ROM:\nα = $(alpha_rom)ppm/°C")
    println(io, "Turner:\nα = $(alpha_turner)ppm/°C")
    println(io, "Kerner:\nα = $(alpha_kerner)ppm/°C")
    println(io, "Shapery:\nα₁ = $(alpha1_shapery)ppm/°C\nα₂ = $(alpha2_shapery)ppm/°C")
    println(io, "Chow:\nα₁ = $(alpha1_chow)ppm/°C\nα₂ = $(alpha2_chow)ppm/°C\nα₃ = $(alpha3_chow)ppm/°C")
    println(io, "Mori Tanaka:\nα₁ = $(alpha1_mt)ppm/°C\nα₂ = $(alpha2_mt)ppm/°C\nα₃ = $(alpha3_mt)ppm/°C")
    
    return String(take!(io))
end


end #module