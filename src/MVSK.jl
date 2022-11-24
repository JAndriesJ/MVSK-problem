module MVSK

include("stock_data.jl")            ; using .stock_data          
include("λ_char.jl")                ; using .λ_char   
include("spaces.jl")                ; using .spaces
include("obj.jl")                   ; using .obj
include("pareto_set.jl")            ; using .pareto_set
include("plot_pareto.jl")           ; using .plot_pareto


# function test_MVSK()
#     include(pwd()*"/test/runtests.jl") 
# end

function get_processed_stock_data()
    return stock_data.read_proc_data()
end

function get_hyperparameter_spaces(mesh_fineness)
    return spaces.get_λ_spaces(mesh_fineness, pd.R_max_max)
end

function get_F_λ_optimal_(M,V,S,K,λ,box_size=0)
    return obj.get_F_λ_opt(M,V,S,K,λ; box_size=box_size, sub=0)
end

function get_F_Λ_Pareto(save_path="default.csv",mesh_fineness=10,box_size=0,sub=0)
    pd = get_processed_stock_data()
    hps =  spaces.get_λ_spaces(mesh_fineness, pd.R_max_max)
    pareto_set.populate_results_csv(pd.R, pd.M, pd.V, pd.S, pd.K, hps.simp_λ_set, hps.simp_conv_mask; save_path=save_path, box_size=box_size, sub=sub)
    return pwd()*"\\"*save_path
end

function plot_Pareto(load_path::String; sel=[1], top_frac=1)
    b         = 0.4433406024717046
    df_pareto = pareto_set.read_results_csv(load_path)
    F         = [df_pareto.F1, df_pareto.F2, df_pareto.F3, df_pareto.F4]
    hps       =  spaces.get_λ_spaces(40, b)
    P         = plot_pareto.plot_Λ(hps, F; sel=sel, top_frac=top_frac)
    return P
end

end


