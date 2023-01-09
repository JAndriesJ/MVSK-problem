module plot_pareto

using PlotlyJS
using Colors, Plots

include("λ_char.jl")                ; using .λ_char 

export plot_Λ,
       plot_moments_vs_λ,
       scale_to_1_optimals

function plot_Λ(hps, f; sel=[1], top_frac =1)
    f = scale_to_1_optimals(f)

    Vol = calc_volume(hps, f, sel, top_frac=top_frac)
    # Sur = calc_surface(hps.tics)
    L = def_Layout(sel)

    P = PlotlyJS.plot(Vol, L)  
    # P = PlotlyJS.plot([Sur, Vol], L)  
    return P
end   

function plot_moments_vs_λ(hps, f)
    moment_names = ["Mean" "Variance"; "Skewness" "Kurtosis"]
    # L = def_Layout([1])
    fig = make_subplots(
        rows=2, cols=2,
        specs=fill(Spec(kind="scene"), 2,2),
        subplot_titles=moment_names,
        vertical_spacing=0.01,
        horizontal_spacing=0.01
    )
   
    for rc in CartesianIndices((2, 2))
        row, col = rc.I
        println([row+col])
        Vol = calc_volume(hps, f, [row+col])
        add_trace!(fig,Vol,row=row, col=col)
    end

    return fig
end  

function calc_volume(hps, f, sel; top_frac =1)
    plot_vals = calc_plot_vals(hps, sel, f)
    isomin = maximum(plot_vals) - maximum(plot_vals)*top_frac
    return  PlotlyJS.volume(x=hps.λ2[:],
                            y=hps.λ4[:],
                            z=hps.λ3[:],
                            value=plot_vals[:],
                            isomin=isomin, #,
                            # isomax=1.0,
                            opacity=0.2, # needs to be small to see through all surfaces
                            surface_count=20, # needs to be a large number for good volume rendering
                            colorscale = "Jet",
                            )
end

function calc_plot_vals(hps, sel, f)
    gd = hps.tics
    plot_vals = -1*ones(gd^3)
    if length(sel) == 1
        plot_vals[hps.simp_mask[:]] = f[sel[1]] 
    else
        plot_vals[hps.simp_mask[:]] = sum(f[sel[i]] for i in 1:length(sel))
    end
    plot_vals = reshape(plot_vals, gd, gd, gd)
end 

function calc_surface(gd)
    data = range(0, stop=1, length=gd)
    z_data = [ ((λ₂+λ₄) ≤ 1.0) && (λ₂+λ₄+λ_char.λ₃_bound_Δ(λ₂, λ₄, 0.4433406)  ≤ 1.05) ? λ_char.λ₃_bound_Δ(λ₂, λ₄, 0.4433406) : NaN for λ₂ in data, λ₄ in data]
    const_color = Plots.cgrad( [ RGB{Float64}(0.5,1,0.5) for _ in 1:2 ] )
    return PlotlyJS.surface(;z=z_data, x=data, y=data, color=const_color, opacity=0.4)
end 

function def_Layout(sel)
    return PlotlyJS.Layout(;
                            title="Showing the values of f$(sel)",
                            xaxis_title = "lambda_2",
                            yaxis_title = "lambda_4",
                            zaxis_title = "lambda_3",
                            scene_camera_eye=attr(x=1, y=-2.5, z=0.25),
                            font=attr(  family="Courier New, monospace",
                                        size=18,
                                        color="RebeccaPurple"
                                    )
                           )

end

function scale_to_1_optimals(arr, flip = false)
    mm = minimum(filter(!isnan, arr))
    mx = maximum(filter(!isnan, arr))
    scal = (arr .- mm)/(mx - mm)
    return flip ? 1 .+ (-1)*scal : scal
end
scale_to_1_optimals(f) = [scale_to_1_optimals(f[1], false),
                          scale_to_1_optimals(f[2], true),
                          scale_to_1_optimals(f[3], false),
                          scale_to_1_optimals(f[4], true)]


# plot_supp_hist(df) =  PlotlyJS.plot(df, x=:supp_w, kind="histogram")

end