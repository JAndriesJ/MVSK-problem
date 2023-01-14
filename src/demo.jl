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
λ = spaces.gen_simp_ele(4)
obj.get_F_λ_opt(pd,λ,x₀=zeros(20),iter=0, box_size=0, sub=0)


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


## Computations
pd  = MVSK.get_processed_stock_data()
hps =  spaces.get_λ_spaces(40, pd.R_Δ_up, pd.R_□_up)

λ_set = hps.simp_λ_set[1:10]
conv_mask = hps.simp_conv_mask[1:10]

pareto_set.prep_results_csv(pd, λ_set, conv_mask,
                            "FISTA_simp_40.csv",
                            iter=2000 ,box_size=0, sub=0)

pareto_set.prep_results_csv(pd, λ_set, conv_mask,
                            "FISTA_box_40.csv",
                            iter=2000 ,box_size=1, sub=0)

pareto_set.prep_results_csv(pd, λ_set, conv_mask,
                            "IPopt_test.csv",
                            iter=0 ,box_size=0, sub=0)
####

df_FISTA = pareto_set.read_results_csv("FISTA_test.csv","|")
df_IPOPT = pareto_set.read_results_csv("IPopt_test.csv","|")

hcat(df_FISTA.Fλ, df_IPOPT.Fλ)
sum(df_FISTA.Fλ .≥ df_IPOPT.Fλ) 
hcat(df_FISTA.t, df_IPOPT.t)
using LinearAlgebra
norm.(df_FISTA.w - df_IPOPT.w)
# .≥ df_IPOPT.w

df_FISTA.t'
df_IPOPT.t'




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
load_path = "C:\\Users\\jandr\\code_projects\\MVSK.jl\\assets\\"* ["pareto_071122_40_simp_sparse_5.csv",
                                                                "pareto_261022_40_simp.csv",
                                                                "pareto_271022_40_box.csv",
                                                                "simp19.12.2022.csv",
                                                                "Box19.12.2022.csv"][4]


df_pareto = pareto_set.read_results_csv(load_path,"\t")

lambdas_df = df_pareto[:,:λ]
w_lambdas_df = df_pareto.w

lambdas = lambdas_df[1:1000]
w_lambdas = w_lambdas_df[1:1000]
ord = sortperm(lambdas)
lambdas = lambdas[ord]
w_lambdas = w_lambdas[ord]


maximum(norm.(lambdas[1:end-1] - lambdas[2:end]))
minimum(norm.(lambdas[1:end-1] - lambdas[2:end]))

using LinearAlgebra
λ_dist_mat = [ norm(l1 - l2, 2) for l1 in lambdas, l2 in lambdas]
wλ_dist_mat = [ norm(w1 - w2, 2) for w1 in w_lambdas , w2 in w_lambdas ]
dily = wλ_dist_mat .≤  8*λ_dist_mat
nara = map(x-> x.I, findall(.~ dily))


using PlotlyJS
plot(heatmap(z=λ_dist_mat))
plot(heatmap(z=wλ_dist_mat))
plot(heatmap(z=(dily.+ 0)))


