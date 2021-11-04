module MVSK

using Pkg
Pkg.status()

include("src\\SDPmodel.jl")
include("src\\SDPoptimized.jl")
using .SDPmodel 
using .SDPoptimized


## Running a single model:
Number_of_stocks = 5
Level_of_Hierarchy = 2 # must be higher than 1
k = 3 # 3 for skewnes, 4 for kurtosis
mod     = SDPmodel.get_SDP_model(Number_of_stocks,Level_of_Hierarchy,3)
mod_opt = SDPoptimized.optimize_SDP(mod)


using JuMP
Primal_status    = string(JuMP.primal_status(mod_opt)) 
Dual_status      = string(JuMP.dual_status(mod_opt))
objective_value  = JuMP.objective_value(mod_opt)
computation_time = JuMP.solve_time(mod_opt)


## Running the model for a batch of parameters:
## Skewness
Number_of_stocks_list = [5:15 ...]
Level_of_Hierarchy_list = [2:3 ...] # must be higher than 1
k = 3 # 3 for skewnes, 4 for kurtosis
SDPoptimized.batch_optimize_SDP(Number_of_stocks_list, Level_of_Hierarchy_list, k)
## Kurtosis
Number_of_stocks_list = [5:15 ...]
Level_of_Hierarchy_list = [2:3 ...] # must be higher than 1
k = 4 # 3 for skewnes, 4 for kurtosis
SDPoptimized.batch_optimize_SDP(Number_of_stocks_list, Level_of_Hierarchy_list, k)




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
