module Model

include("stock_data.jl")   ; using .stock_data 
include("data_moments.jl") ; using .data_moments
include("convex_quadratic_optimization.jl")  ; using .convex_quadratic_optimization
#include("objective_functions.jl") ; using .objective_functions

export build_model,
       test


function build_model(num_stocks = 6, λ = ones(4)*0.25 )
    R = stock_data.load_cont_comp_returns()[:,1:num_stocks]
    moment_matrices = Dict( "means_vector"      => data_moments.get_means_vector(R), 
                            "covariance_matrix" => data_moments.get_covariance_matrix(R),
                            "skewness_matrix"   => data_moments.get_skewness_matrix(R),
                            "kurtosis_matrix"   => data_moments.get_kurtosis_matrix(R))
    f₁ₘₐₓ = maximum(moment_matrices["means_vector"])

    f₂ₘᵢₙ = convex_quadratic_optimization.get_variance_min(moment_matrices["covariance_matrix"]*(10^4))*(10^(-4)) # scaling for stability
    return MVSK_data(num_stocks, R , moment_matrices, f₁ₘₐₓ, f₂ₘᵢₙ, λ)
end

struct MVSK_data
    n_stocks::Int64
    Data_matrix::Matrix{Float64}
    moment_matrices::Dict{String, Matrix{Float64}}
    f₁ₘₐₓ::Float64
    f₂ₘᵢₙ::Float64
    λ::Vector{Float64}
    #stock_names::Vector{String}
    # begin_date::Dates.Date
    # end_date::Dates.Date
end

function test()
    data_moments.run_tests()
    convex_quadratic_optimization.run_tests()

    error("There are no tests for stock data")
end

end