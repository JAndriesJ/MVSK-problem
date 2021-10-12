module MVSK

# using Pkg
# Pkg.status()

include("src\\stock_data.jl")
include("src\\moments.jl")
include("src\\constraints.jl")
include("src\\SDPmodel.jl")
include("src\\SDPoptimized.jl")
using .stock_data
using .moments 
using .constraints 
using .SDPmodel 
using .SDPoptimized




function batchrun()
    t = 3 
    for N in 10:25 
        SDPmod = SDPmodel.get_SDP_model(N,t,4)
        optSDPmod = SDPoptimized.optimize_SDP(SDPmod)
        abc = open("C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\Computation times.txt", "a") 
            write(abc, "Computed kurtosis maximization bound: $(JuMP.objective_value(optSDPmod)) at level $(t) for $(N)-stockstook $(JuMP.solve_time(optSDPmod)) seconds"*"\n") 
        close(abc)
    end
end

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
