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


using JuMP


t = 3 
for N in 10:25 
    SDPmod = SDPmodel.get_SDP_model(N,t,4)
    optSDPmod = SDPoptimized.optimize_SDP(SDPmod)
    abc = open("C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\Computation times.txt", "a") 
        write(abc, "Computed kurtosis maximization bound: $(JuMP.objective_value(optSDPmod)) at level $(t) for $(N)-stockstook $(JuMP.solve_time(optSDPmod)) seconds"*"\n") 
    close(abc)
end



# line_by_line = readlines("C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\S&P500Tickers.csv")

# abc = open("C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\123.csv", "a") 
# for line in line_by_line 
#     write(abc, replace(line[2:end-4], "\"" => "")*"\n") 
# end
# close(abc)


mutable struct Bound_data{N,lvl,k}
    Number_of_stocks::N
    level::lvl
    objective::k
end
bound_data10_2_3 = Bound_data{10,2,3}


tickers = stock_data.get_tickers()
load_relative_returns_data_matrix()
R = load_relative_returns_data_matrix([1:30...])
α = [1,0,0,0,0,2,0,0,0,0]
get_data_moment(R,α)


return
end

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
