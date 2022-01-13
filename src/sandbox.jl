## Initialize 
using Pkg
Pkg.activate(dirname(@__DIR__))

include("MVSK.jl") ; using .MVSK

MVSK_data = MVSK.get_MVSK_data()
MVSK_bounds = MVSK.get_MVSK_bounds(MVSK_data)

print(fieldnames(typeof(MVSK_data)))
MVSK_data.n_stocks
MVSK_data.Data_matrix
MVSK_data.moment_matrices


print(fieldnames(typeof(MVSK_bounds)))
MVSK_bounds.λ
MVSK_bounds.deg_range
MVSK_bounds.grid_range

MVSK_bounds.f₁ₘₐₓ
MVSK_bounds.f₂ₘᵢₙ
MVSK_bounds.f₃ₘₐₓ
MVSK_bounds.f₄ₘᵢₙ
MVSK_bounds.𝒻_tilde

MVSK.save_MVSK(MVSK_data, MVSK_bounds)



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


