
using Pkg
Pkg.status()
Pkg.activate("C:\\Users\\jandr\\MahCodes\\MVSK\\")

include("stock_data.jl")
include("Data_moments.jl")
include("CQP_simplex.jl")
include("Grid_Search.jl")
include("POP\\using_SumOfSquares.jl")
using .stock_data 
using .Data_moments
using .CQP_simplex
using .Grid_Search 
using .using_SumOfSquares

Data_moments.run_tests()
CQP_simplex.run_tests()
Grid_Search.run_tests()
using_SumOfSquares.run_tests()

n = 6
R = stock_data.load_relative_centralize_returns_data_matrix()[:,1:n]

## Estimating f₁ₘₐₓ and f₂ₘᵢₙ
f₁ₘₐₓ =  maximum(sum(R, dims= 1))/(size(R)[2]-1)
covariance_matrix = Data_moments.get_covariance_matrix(R)
f₂ₘᵢₙ = CQP_simplex.get_variance_min(covariance_matrix)  

## Estimating f₃ₘₐₓ and f₄ₘᵢₙ 
## Towards lower bounds of 
f₃  = Data_moments.get_f₃(R)
f₄  = Data_moments.get_f₄(R)
# via Grid search
r = 4
hat_f₃ₘₐₓ = Grid_Search.max(f₃,n,r)
hat_f₄ₘᵢₙ = Grid_Search.min(f₄,n,r)
# via interior point line search filter method

# via multi multi-start
Data_moments.get_means_vector(R*10^6)
# via particle swarm

## Towards upper bounds of 
# via SumsOfSquares.jl
skewness_matrix = Data_moments.get_skewness_matrix(R)
tilde_f₃ₘₐₓ = -using_SumOfSquares.get_SOS_bound((-1)*skewness_matrix,3)
kurtosiss_matrix = Data_moments.get_kurtosis_matrix(R)
tilde_f₄ₘᵢₙ = maximum([using_SumOfSquares.get_SOS_bound(kurtosiss_matrix,6),1e-16]) # WARNING: somewhat arbitary value
using_SumOfSquares.get_SOS_bound(kurtosiss_matrix,4)
using_SumOfSquares.get_SOS_bound(kurtosiss_matrix,6)
# via (my coded) Lasser Hierarchy (there are some complications)
# include("POP\\Lasserre_hier\\Lasserre_bound.jl")
# using .Lasserre_bound 
# Lasserre_bound.get_Lasserre_bound(-skewness_matrix,3)


## Bounding funcitons for 𝒻:
λ = [0.25,0.25,0.25,0.25]
### Upper bounding function
tilde_ϕ = f₁ₘₐₓ,f₂ₘᵢₙ,tilde_f₃ₘₐₓ,tilde_f₄ₘᵢₙ
tilde_𝒻 = Data_moments.get_𝒻(λ,tilde_ϕ,R)

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
