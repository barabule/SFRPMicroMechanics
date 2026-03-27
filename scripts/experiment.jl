

module Experiment

using GLMakie
import Makie: Relative, Auto, Consume

function LogRelativeSlider(fig_or_pos; 
                startvalue = 1.0, 
                range = (1e-3, 1e6), 
                label_prefix = "Value: ",
                width = 300.0,
                height = 30.0,
                bg_color = :gray95,
                strokecolor = :gray80,
                drag_active_color = :skyblue,
                drag_deactive_color = :lightgray,
                )
    # Create a nested layout in the target slot to manage the stacking of Box/Label/Textbox
    gl = fig_or_pos[] = GridLayout()
    
    # State Observables
    val = Observable(Float64(startvalue))
    dragging = Observable(false)
    editing = Observable(false)
    fill_p = Observable(0.5)
    
    # 1. Background Track
    bg = Box(gl[1, 1], color = bg_color, strokecolor = strokecolor)
    
    # 2. The Fill Bar (Visible only when not editing)
    bar = Box(gl[1, 1], 
        color = lift(d -> d ? drag_active_color : drag_deactive_color, dragging),
        width = lift(p -> Relative(p), fill_p),
        halign = :left,
        visible = lift(!, editing)
    )
    
    # 3. The Text Label (Visible only when not editing)
    lbl = Label(gl[1, 1], lift(v -> "$label_prefix$(round(v, digits=3))", val),
            visible = lift(!, editing),
            halign = :left,        # Align text to the left of the box
            padding = (30, 0, 0, 0) # Give it a little breathing room from the edge
        )
    
    # WORKAROUND: Explicitly use RGBAf for all color states to avoid MethodErrors
    tbox = Textbox(gl[1, 1], 
        placeholder = "Type...",
        width = lift(e -> e ? width : 0.0, editing), 
        height = lift(e -> e ? height : 0.0, editing),
        textcolor = lift(e -> e ? RGBAf(0, 0, 0, 1) : RGBAf(0, 0, 0, 0), editing),
        boxcolor = lift(e -> e ? RGBAf(1, 1, 1, 1) : RGBAf(0, 0, 0, 0), editing),
        bordercolor = lift(e -> e ? RGBAf(0.5, 0.5, 0.5, 1) : RGBAf(0, 0, 0, 0), editing),
        halign = :center,
        valign = :center
    )

    on(tbox.stored_string) do s
        parsed = tryparse(Float64, s)
        if !isnothing(parsed)
            val[] = clamp(parsed, range[1], range[2])
        end
        editing[] = false
        try tbox.focused[] = false catch; end # Release focus
    end

    # Event Interaction Logic
    parent_scene = Makie.get_topscene(fig_or_pos)
    last_click_time = Ref(0.0)
    last_mouse_pos = Ref(Point2f(0))
    
    on(events(parent_scene).mousebutton) do event
        # Detect if mouse is inside the widget using 'computedbbox'
        bbox = bg.layoutobservables.computedbbox[]
        mouse_pos = events(parent_scene).mouseposition[]
        mouse_inside = mouse_pos ∈ bbox

        if event.button == Mouse.left
            if event.action == Mouse.press && mouse_inside
                t = time()
                # Double Click Detection
                if t - last_click_time[] < 0.3
                    editing[] = true
                    tbox.displayed_string[] = string(round(val[], digits=4))
                    # Trigger focus so the user can type immediately
                    try tbox.focused[] = true catch; end 
                else
                    if !editing[]
                        dragging[] = true
                        last_mouse_pos[] = mouse_pos
                    end
                end
                last_click_time[] = t
            elseif event.action == Mouse.release
                dragging[] = false
                fill_p[] = 0.5 # Visual "joystick" reset
            end
        end
        # Consume events if we are typing so we don't trigger other shortcuts
        return Consume(editing[])
    end
    on(events(parent_scene).mouseposition) do pos
    if dragging[] && !editing[]
        delta = pos[1] - last_mouse_pos[][1]
        sensitivity = ispressed(parent_scene, Keyboard.left_shift) ? 0.0005 : 0.005
        
        # 1. Get current value
        curr = val[]
        
        # 2. Hybrid Logic: 
        # If the value is tiny or zero, we use a small linear step to "kickstart" it.
        # Otherwise, we use the logarithmic scaling.
        epsilon = 1e-3 
        
        if abs(curr) < epsilon
            # Linear "kick" to get away from zero
            new_val = curr + (delta * sensitivity)
        else
            # Standard Logarithmic scaling
            # We use sign(curr) so it works for negative ranges too!
            new_val = curr * 10^(delta * sensitivity * sign(curr))
        end
        
        val[] = clamp(new_val, range[1], range[2])
        fill_p[] = clamp(fill_p[] + (delta / 400.0), 0.0, 1.0)
        last_mouse_pos[] = pos
    end
    return Consume(false)
end

    on(events(parent_scene).mouseposition) do pos
        if dragging[] && !editing[]
            delta = pos[1] - last_mouse_pos[][1]
            
            # Sensitivity multipliers: Shift for fine-tuning
            sensitivity = ispressed(parent_scene, Keyboard.left_shift) ? 0.0005 : 0.005
            
            # Apply Logarithmic Update
            new_val = val[] * 10^(delta * sensitivity)
            val[] = clamp(new_val, range[1], range[2])
            
            # Update Visual Fill
            fill_p[] = clamp(fill_p[] + (delta / 400.0), 0.0, 1.0)
            last_mouse_pos[] = pos
        end
        return Consume(false)
    end

    return (value = val, layout = gl)
end

end #module



function main()
fig = Experiment.Figure()
# A slider that stays between 0.1 and 100.0
s = Experiment.LogRelativeSlider(fig[1, 1], startvalue = 1.0, range = (0.1, 100.0))

# Display the value elsewhere just to verify
Experiment.Label(fig[2, 1], Experiment.GLMakie.lift(v -> "External Listener: $v", s.value))

fig

end



function main2()

    fig = Experiment.Figure(size = (400, 300))

    # 1. Define your parameters
    params = [
        (name = "Frequency", start = 10.0,  range = (1.0, 100.0)),
        (name = "Amplitude", start = 1.0,   range = (0.1, 10.0)),
        (name = "Phase Shift", start = 0.01, range = (0.001, 1.0)),
        (name = "Decay Rate", start = 0.5,  range = (0.01, 2.0))
    ]

    # 2. Create a container for the sliders
    # Using a nested layout makes it easy to move the whole block later
    controls_layout = fig[1, 1] = Experiment.GridLayout()

    # Store the observables so you can use them in your plots
    outputs = Dict()

    for (i, p) in enumerate(params)
        # Create the widget in row 'i'
        # We use label_prefix to show the name
        slider_unit = Experiment.LogRelativeSlider(
            controls_layout[i, 1], 
            startvalue = p.start, 
            range = p.range, 
            label_prefix = "$(p.name): "
        )
        
        outputs[p.name] = slider_unit.value
    end

    # 3. UI Polish
    Experiment.rowgap!(controls_layout, 15) # Tighten the vertical space between sliders

    # If you want the sliders to have a specific height
    # rowsize!(controls_layout, 1, Fixed(30)) 

    fig

end


function main3()


    function decaying_sine(t, p)
        f, a, ϕ, d = p
        return a .* sin.(2π .* f .* t .+ ϕ) .* exp.(-d .* t)
    end

    # -- Setup the Main Figure --
    fig = Experiment.Figure(size = (1000, 600))
    Experiment.colsize!(fig.layout, 1, Experiment.Relative(0.3)) # Allocate 30% width for controls

    # -- 1. The Sliders (Controls) --
    controls_grid = fig[1, 1] = Experiment.GridLayout( tellheight = false)

    # Define parameters for 4 different controls
    params = [
        (name = "Frequency", start = 1.0,  range = (0.1, 10.0)),
        (name = "Amplitude", start = 5.0,  range = (0.1, 10.0)),
        (name = "Phase Shift", start = 0.0, range = (-1.0, 1.0)),
        (name = "Decay Rate", start = 0.2,  range = (0.01, 2.0))
    ]

    # We need a Dict to store the value Observables from our function
    outputs = Dict()

    for (i, p) in enumerate(params)
        # Create the widget and place it in a new row (i)
        slider_data = Experiment.LogRelativeSlider(
            controls_grid[i, 1], 
            startvalue = p.start, 
            range = p.range, 
            label_prefix = "$(p.name): "
        )
        # Store the .value Observable
        outputs[p.name] = slider_data.value
    end

    # Fine-tune the control grid layout
    Experiment.rowgap!(controls_grid, 5) # No vertical gap for a tight list
    # Experiment.colsize!(controls_grid, 1, Experiment.Fixed(800))
    # -- 2. The Master Observable and Plot --
    ax = Experiment.Axis(fig[1, 2], title = "Decaying Sine Control", xlabel = "Time", ylabel = "Amplitude")
    # Set a fixed range for a better "Blender-like" feeling
    Experiment.ylims!(ax, -10, 10) 

    # Create the time vector
    t = range(0, 10, length=1000)

    # The MAGIC: Combine all 4 observables into a single "state" tuple
    all_values = Experiment.lift(outputs["Frequency"], outputs["Amplitude"], outputs["Phase Shift"], outputs["Decay Rate"]) do f, a, p, d
        return (f, a, p, d) # This tuple updates whenever ANY value changes
    end

    # Lift the plotting function to update whenever the state changes
    data = Experiment.lift(p -> decaying_sine(t, p), all_values)

    # Add the line plot
    Experiment.lines!(ax, t, data, color = :royalblue, linewidth = 2)

    # Set the final layout
    Experiment.colgap!(fig.layout, 10) # Simple gap between controls and plot

    fig

end