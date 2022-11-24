include("MVSK.jl")            ; using .MVSK     

#-------------------------------------------- Data -------------------------------------------- 
pd = MVSK.get_processed_stock_data()
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


#-------------------------------------------- -------------------------------------------- 
# MVSK.get_F_Λ_Pareto(save_path="default.csv",mesh_fineness=10,box_size=0,sub=0)
load_path = "C:\\Users\\jandr\\code_projects\\MVSK\\assets\\"* ["pareto_071122_40_simp_sparse_5.csv",
                                                                "pareto_261022_40_simp.csv",
                                                                "pareto_271022_40_box.csv"][2]
MVSK.plot_Pareto(load_path, sel=[4])


