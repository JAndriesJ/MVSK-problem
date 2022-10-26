module MVSK

include("stock_data.jl")            ; using .stock_data          
include("λ_char.jl")                ; using .λ_char   
include("spaces.jl")                ; using .spaces
include("obj.jl")                   ; using .obj
include("pareto_set.jl")            ; using .pareto_set
include("plot_pareto.jl")           ; using .plot_pareto

# include(pwd()*"/test/runtests.jl")  
#----------------------------------------------------------------------------------------Data 
pd  = stock_data.get_proc_data()
hps =  spaces.get_λ_spaces(10, pd.R_max_max)
save_path = "C:\\Users\\jandr\\code_projects\\MVSK\\assets\\pareto_1.csv"
pareto_set.populate_results_csv(pd.R, pd.M, pd.V, pd.S, pd.K, hps.simp_λ_set, hps.simp_conv_mask; save_path = save_path)

df_pareto = pareto_set.read_results_csv(save_path)
f = [df_pareto.f1, df_pareto.f2, df_pareto.f3, df_pareto.f4]




f[1][calc_close_to_best(f[1])]
f[2][calc_close_to_best(f[2], flip = true)]
f[3][calc_close_to_best(f[3])]
f[4][calc_close_to_best(f[4], flip = true)]

nara = scale_to_1_optimals(f[1])+ scale_to_1_optimals(f[2],true)+ scale_to_1_optimals(f[3]) + scale_to_1_optimals(f[4],true)
best_mast = nara .>  maximum(nara)-0.001*maximum(nara)

best_F_mat = hcat(f[1],f[2],f[3],f[4])[best_mast,:]




function scale_to_1_optimals(arr, flip = false)
    mm = minimum(filter(!isnan, arr))
    mx = maximum(filter(!isnan, arr))
    scal = (arr .- mm)/(mx - mm)
    return flip ? 1 .+ (-1)*scal : scal
end



using PlotlyJS
P = plot_pareto.plot_Λ(hps, f; sel=[1,2,3])
PlotlyJS.savefig(P, "test.png")



# https://plotly.com/julia/figure-labels/
# https://plotly.com/julia/3d-surface-plots/
end


# sort(f[1])[end-4:end]'
# sort(f[2])[1:5]'
# sort(f[3])[end-4:end]'
# sort(f[4])[1:5]'

function calc_close_to_best(arr; closenes = 0.1, flip = false)
    s_arr = scale_to_1_optimals(arr, flip)
    return s_arr .>  1-closenes
end
