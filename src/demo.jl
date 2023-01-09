# using Pkg
# Pkg.status()
# Pkg.activate(normpath(joinpath(@__DIR__, "..")))
src_p = ""
include("MVSK.jl")                        ; using .MVSK     
include(src_p*"stock_data.jl")            ; using .stock_data          
include(src_p*"λ_char.jl")                ; using .λ_char   
include(src_p*"spaces.jl")                ; using .spaces
include(src_p*"obj.jl")                   ; using .obj
include(src_p*"pareto_set.jl")            ; using .pareto_set
include(src_p*"plot_pareto.jl")           ; using .plot_pareto

#-------------------------------------------- Data -------------------------------------------- 
pd = MVSK.get_processed_stock_data()
hps = MVSK.get_hyperparameter_spaces(80, pd.R_Δ_up, pd.R_□_up)
approximate_volumns = round.(hps.conv_rel_vol, digits=3)

### Plots
using PlotlyJS
using Colors

plot_vals = hps.simp_mask .+ (hps.simp_mask .* hps.conv_plus_mask) .+  (hps.simp_mask .* hps.conv_□_mask) .+  (hps.simp_mask .* hps.conv_Δ_mask)
Vol =  PlotlyJS.volume( x=hps.λ2[:],
                        y=hps.λ4[:],
                        z=hps.λ3[:],
                        value=plot_vals[:],
                        isomin=0.1, #,
                        # isomax=0.8,
                        opacity=0.2, # needs to be small to see through all surfaces
                        surface_count=17, # needs to be a large number for good volume rendering
                        colorscale = "Jet"
                        )
PlotlyJS.plot(Vol)



λ_ =  range(0, stop=1, length=100)
λ_grid = [[λ4,λ3] for λ4 in λ_, λ3 in λ_ ]
nar = [ sum(λ_g) ≤ 1  for λ_g ∈ λ_grid]
nar1 = [ sum(λ_g) ≤ 1 ? λ_g[2] ≤ sqrt(abs((8/3)*λ_g[1]*(1-λ_g[1]-λ_g[2]))) : 0   for λ_g ∈ λ_grid]
# nar2 = [ λ_g[2] ≤ sqrt((8/3)(1-λ_g[1]-λ_g[2] )*λ_g[1] )  for λ_g ∈ λ_grid]
# pd.R_Δ_up, pd.R_□_up


using PlotlyJS

data = [1 20 30; 20 1 60; 30 60 1]
plot(heatmap(z=nar*.1 + nar1*.1 ))



## Computations
# MVSK.get_F_Λ_Pareto("simp19.12.2022.csv",40,0,0)
# MVSK.get_F_Λ_Pareto("Box19.12.2022.csv",40,1,0)
pd        = MVSK.get_processed_stock_data()
hps =  spaces.get_λ_spaces(40, pd.R_Δ_up, pd.R_□_up)
pareto_set.populate_results_csv(pd.R, pd.M, pd.V, pd.S, pd.K, hps.simp_λ_set[10567:end], hps.simp_conv_mask[10567:end]; save_path="Box19.12.2022.csv", box_size=1, sub=0)


## Post proc
load_path = [ "simp19.12.2022.csv", "Box19.12.2022.csv"][1]
df_pareto = pareto_set.read_results_csv(load_path,"\t")
F         = [df_pareto.F1, df_pareto.F2, df_pareto.F3, df_pareto.F4]
F_sca     = scale_to_1_optimals(F)
F_sca_sum = F_sca[1] .+ F_sca[2] .+ F_sca[3] .+ F_sca[4] 
max_F_sca_sum = maximum(F_sca_sum)
top_1_mask =  F_sca_sum .> 0.99*max_F_sca_sum
F_top_1 = [F[1][top_1_mask], F[2][top_1_mask], F[3][top_1_mask], F[4][top_1_mask]]
F_top_1[1]
F_top_1[4]

minimum(F[3]),maximum(F[3])
minimum(F_top_1[3]),maximum(F_top_1[3])


hps       =  spaces.get_λ_spaces(40, pd.R_Δ_up, pd.R_□_up)

plot_pareto.plot_Λ(hps, F; sel=[1], top_frac=1)
plot_pareto.plot_Λ(hps, F; sel=[2], top_frac=1)
plot_pareto.plot_Λ(hps, F; sel=[3], top_frac=1)
plot_pareto.plot_Λ(hps, F; sel=[4], top_frac=1)

plot_pareto.plot_Λ(hps, F; sel=[1,2,3,4], top_frac=0.025)





#--------------------------------------------single calc-------------------------------------------- 
λ = [0.2878  0.302  0.2831  0.1271]
simp_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,λ)
simp_F_λ_opt.domain
simp_F_λ_opt.optimizer 
simp_F_λ_opt.opt_val
simp_F_λ_opt.opt_stat
simp_F_λ_opt.sol_time
simp_F_λ_opt.λ

box_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,λ,1)
box_F_λ_opt.domain
box_F_λ_opt.optimizer 
box_F_λ_opt.opt_val
box_F_λ_opt.opt_stat
box_F_λ_opt.sol_time
box_F_λ_opt.λ

# Sparse one ??
simp_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,Λ[1])

a = zeros(10)
Λ = [ spaces.gen_simp_ele() for i ∈ 1:10] 

Threads.@threads for i = 1:2
    λ = Λ[i]
    # simp_F_λ_opt = MVSK.get_F_λ_optimal_(pd.M,pd.V,pd.S,pd.K,λ)
    # a[i] = simp_F_λ_opt.sol_time
    a[i] = Λ[i]
end


#-------------------------------------------- -------------------------------------------- 
# MVSK.get_F_Λ_Pareto(save_path="default.csv",mesh_fineness=10,box_size=0,sub=0)
load_path = "C:\\Users\\jandr\\code_projects\\MVSK\\assets\\"* ["pareto_071122_40_simp_sparse_5.csv",
                                                                "pareto_261022_40_simp.csv",
                                                                "pareto_271022_40_box.csv"][2]
MVSK.plot_Pareto(load_path, sel=[4])


f₁(w )
f₂(w )
f₃(w )
f₄(w )


a = zeros(10)
Threads.@threads for i = 1:10
    a[i] = Threads.threadid()
end
a
Threads.nthreads()
