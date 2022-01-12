## Initialize 
using Pkg
Pkg.status()
Pkg.activate(dirname(@__DIR__))

############################ MODEL ########################################
include("MODEL\\Model.jl")   ; using .Model 
mod = build_model(6)

# mod.n_stocks
# mod.Data_matrix
# mod.moment_matrices
# mod.f₁ₘₐₓ
# mod.f₂ₘᵢₙ 

############################## POP ######################################
me_vec  = mod.moment_matrices["means_vector"]
cov_mat = mod.moment_matrices["covariance_matrix"]
sk_mat  = mod.moment_matrices["skewness_matrix"]
kur_mat = mod.moment_matrices["kurtosis_matrix"]


include("POP\\SOS.jl") ; using .SOS
deg_range = 4:2:8
# https://jump.dev/SumOfSquares.jl/stable/generated/Polynomial%20Optimization/polynomial_optimization/#The-maxdegree-keyword-argument

f₃ₘₐₓ_pop_df = SOS.get_SOS_skewness_estimates(sk_mat, deg_range)
f₄ₘᵢₙ_pop_df = SOS.get_SOS_kurtosis_estimates(kur_mat, deg_range) 

include("MODEL\\convex_quadratic_optimization.jl")  ; using .convex_quadratic_optimization
f₄ₘᵢₙ_cqr = convex_quadratic_optimization.get_relaxed_kurtosis(kur_mat)

################################ OBJ FUNCT ####################################
include("UTIL\\objective_functions.jl"); using .objective_functions
f₁, f₂, f₃, f₄ = objective_functions.get_f_vec(mod.moment_matrices)
𝒻 = objective_functions.get_𝒻(ones(4), 
                              (mod.f₁ₘₐₓ, mod.f₂ₘᵢₙ, f₃ₘₐₓ_pop_df[end,2], abs(f₄ₘᵢₙ_pop_df[end,2])), # fₒₚₜ 
                              mod.moment_matrices)
                           
################################ NON LIN ####################################
grid_range = 2 .^ (1:6)
include("NONLIN\\grid_search.jl")  ; using .grid_search 
gs_upper_bounds = grid_search.batch_min(𝒻, mod.n_stocks, grid_range)

f₃ₘₐₓ_nl_df = grid_search.batch_max(f₃, mod.n_stocks, grid_range)
f₄ₘᵢₙ_nl_df = grid_search.batch_min(f₄, mod.n_stocks, grid_range)

################################ SAGE ####################################
# Store data and switch to python:
using CSV, DataFrames, Dates

function save_data(mod, f₃ₘₐₓ_nl_df, f₄ₘᵢₙ_nl_df, f₃ₘₐₓ_pop_df, f₄ₘᵢₙ_pop_df, gs_upper_bounds; delim=',')
    save_dir = "assets\\"*replace(string(now()), ":" => "-")
    mkdir(save_dir)
    # save moment matrices
    for mom_mat ∈ ["means_vector", "kurtosis_matrix", "covariance_matrix", "skewness_matrix"]
        CSV.write(save_dir*"\\"*mom_mat*".csv", DataFrame(mod.moment_matrices[mom_mat]),
                                            header = false,
                                            delim=delim)
    end

    CSV.write(save_dir*"\\f₃ₘₐₓ_nl_df.csv", f₃ₘₐₓ_nl_df, header = true, delim=delim)
    CSV.write(save_dir*"\\f₄ₘᵢₙ_nl_df.csv", f₄ₘᵢₙ_nl_df, header = true, delim=delim)

    CSV.write(save_dir*"\\f₃ₘₐₓ_pop_df.csv", f₃ₘₐₓ_pop_df, header = true, delim=delim)
    CSV.write(save_dir*"\\f₄ₘᵢₙ_pop_df.csv", f₄ₘᵢₙ_pop_df, header = true, delim=delim)

    CSV.write(save_dir*"\\gs_upper_bounds.csv", gs_upper_bounds, header = true, delim=delim)
    return save_dir
end

save_dir =save_data(mod,
                    f₃ₘₐₓ_nl_df,
                    f₄ₘᵢₙ_nl_df,
                    f₃ₘₐₓ_pop_df,
                    f₄ₘᵢₙ_pop_df,
                    gs_upper_bounds,
                    delim='|')


nar = DataFrame(Dict( "f₁ₘₐₓ" => mod.f₁ₘₐₓ,  
                      "f₂ₘᵢₙ" => mod.f₂ₘᵢₙ,
                      "f₄ₘᵢₙ_cqr" => f₄ₘᵢₙ_cqr
                    ))


CSV.write(save_dir*"\\parameters.csv", nar, header = true, delim="|")


####################################################################

##### Code Health check

# objective_functions.run_tests()
# grid_search.run_tests()
# using_SumOfSquares.run_tests()
  


# via multi multi-start (TODO)

# via particle swarm (TODO)


## Towards upper bounds of 



# via (my coded) Lasser Hierarchy (there are some complications)
# include("POP\\Lasserre_hier\\Lasserre_bound.jl")
# using .Lasserre_bound 

# Lasserre_bound.get_Lasserre_bound(-skewness_matrix,3)


