module MVSK

# using Pkg
# Pkg.status()

include("src\\stock_data.jl")
include("src\\moments.jl")
include("src\\SDPconstraints.jl")
include("src\\SDPmodel.jl")
include("src\\SDPoptimized.jl")
using .stock_data
using .moments 
using .SDPconstraints 
using .SDPmodel 
using .SDPoptimized


# N,t,k = 5,3,3
# SDP_model = get_SDP_model(N,t,k)

N_list,t_list,k = [5:10 ...], [2:3 ...], 3
SDPoptimized.batch_optimize_SDP(N_list,t_list,k)
N_list,t_list,k = [5:10 ...], [2:3...], 4
SDPoptimized.batch_optimize_SDP(N_list,t_list,k)


using Test, DataFrames
using JuMP




mutable struct Bound_data{N,lvl,k}
    Number_of_stocks::N
    level::lvl
    objective::k
end
bound_data10_2_3 = Bound_data{10,2,3}



mutable struct Portfolio <: Any # portfolio_candidates
    number_of_stocks::Int
    number_of_data_point::Int 
    stock_tickers::Vector{String} 
    relative_returns::Matrix{Float64} 
end
ticker_list = get_tickers()
ticker_list[1:5]

P1 = Portfolio(5,99,["A", "AAL", "AAP", "AAPL", "ABBV"],)
P1.number_of_stocks



end
