## Initialize 
using Pkg
Pkg.status()


Pkg.activate(dirname(@__DIR__))

####################################################################
include("MODEL\\Model.jl")   ; using .Model 
mod = build_model(6)

mod.n_stocks
mod.Data_matrix
mod.moment_matrices
mod.f₁ₘₐₓ
mod.f₂ₘᵢₙ

####################################################################

include("POP\\POP.jl")   ; using .POP

sk_mat = mod.moment_matrices["skewness_matrix"]
kur_mat = mod.moment_matrices["kurtosis_matrix"]

sk_df = POP.get_POP_estimates(sk_mat,3:8)
kur_df = POP.get_POP_estimates(kur_mat,4:8)






df

# obj_val, p_stat, d_stat, sol_t =


####################################################################

include("grid_search.jl")  ; using .grid_search 





## Code Health check

objective_functions.run_tests()
grid_search.run_tests()
using_SumOfSquares.run_tests()
  

## Estimating f₃ₘₐₓ and f₄ₘᵢₙ 
f₃  = get_f₃(skewness_matrix)
f₄  = get_f₄(kurtosis_matrix)
### Towards lower bounds of 

# via Grid search
r = 4
hat_f₃ₘₐₓ = grid_search.max(f₃,n,r)
hat_f₄ₘᵢₙ = grid_search.min(f₄,n,r)
# via interior point line search filter method

# via multi multi-start (TODO)

# via particle swarm (TODO)

## Towards upper bounds of 
# via SumsOfSquares.jl
tilde_f₃ₘₐₓ = -using_SumOfSquares.get_SOS_bound((-1)*skewness_matrix,3)
tilde_f₄ₘᵢₙ = maximum([using_SumOfSquares.get_SOS_bound(kurtosiss_matrix,6),1e-16]) # WARNING: somewhat arbitary value
# via (my coded) Lasser Hierarchy (there are some complications)
# include("POP\\Lasserre_hier\\Lasserre_bound.jl")
# using .Lasserre_bound 
# Lasserre_bound.get_Lasserre_bound(-skewness_matrix,3)


## Bounding funcitons for 𝒻:
λ = [0.25,0.25,0.25,0.25]
f_opts = f₁ₘₐₓ, f₂ₘᵢₙ, tilde_f₃ₘₐₓ, tilde_f₄ₘᵢₙ
tilde_𝒻 = objective_functions.get_𝒻(λ, f_opts, R)

a = [1,1,1,1,1,1]/6
b = [0,0,0,0,0,0]
tilde_𝒻(a)

ϕ = f₁ₘₐₓ,f₂ₘᵢₙ,tilde_f₃ₘₐₓ,tilde_f₄ₘᵢₙ
f₁  = Data_moments.get_f₁(R); f₂  = Data_moments.get_f₂(R); f₃  = Data_moments.get_f₃(R); f₄  = Data_moments.get_f₄(R) 
function 𝒻(x)
    if f₃(x) ≤ ϕ[3] && f₄(x) ≥ ϕ[4] 
        return (1 - (f₁(x)/ϕ[1]))^λ[1] + ((f₂(x)/ϕ[2]) - 1)^λ[2] + (1 - (f₃(x)/ϕ[3]))^λ[3] + ((f₄(x)/ϕ[4]) - 1)^λ[4] 
    else 
        return 0
    end
end

### Lower bounding function
hat_ϕ = f₁ₘₐₓ,f₂ₘᵢₙ,tilde_f₃ₘₐₓ,tilde_f₄ₘᵢₙ
hat_𝒻 = Data_moments.get_𝒻(λ,hat_ϕ,R)


𝒻(a)
