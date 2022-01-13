module stock_data
using CSV, DataFrames
data_dir = string(@__DIR__)*"\\assets\\"

export  read_data,
        load_cont_comp_returns,
        load_relative_returns,
        run_tests

global data_dir = pwd()*"\\assets\\"

function read_data()
    csv_file = data_dir*"stock_prices.csv"
    df = CSV.read(csv_file, DataFrame)
    return select!(df[end-500:end,:], Not(:date))
end


#load_relative_centralize_returns_data_matrix() = centralize_returns(load_relative_returns_data_matrix())
load_div_returns() = get_div_returns(Matrix(read_data()))
load_cont_comp_returns() = get_cont_comp_returns(Matrix(read_data()))
load_relative_returns() = get_relative_returns(Matrix(read_data()))

# transformations
get_div_returns(R) = R[2:end,:] ./ R[1:end-1,:]
get_cont_comp_returns(R)  = log.(get_div_returns(R))

get_diff_returns(R) = diff(R,dims=1) 
get_relative_returns(R) = get_diff_returns(R)./ R[1:end-1,:]


function run_tests()
    error("No tests present")
end


end



# load_relative_returns_data_matrix(ticker) = get_relative_returns(Matrix(read_data(ticker)))