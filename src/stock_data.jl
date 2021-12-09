module stock_data
using CSV, DataFrames
data_dir = string(@__DIR__)*"\\assets\\"

export  read_data,
        load_relative_returns_data_matrix,
        load_relative_centralize_returns_data_matrix

global data_dir = pwd()*"\\assets\\"

function read_data()
    csv_file = data_dir*"stock_prices.csv"
    df = CSV.read(csv_file, DataFrame)
    return select!(df[end-500:end,:], Not(:date))
end

read_data(ticker::String) = read_data()[:,Symbol(ticker)]
read_data(ticker::Int) = read_data()[:,ticker]
read_data(ticker_list::Vector{String}) = read_data()[:,Symbol.(ticker_list)]
read_data(ticker_list::Vector{Int}) = read_data()[:,ticker_list]

load_relative_returns_data_matrix(ticker) = get_relative_returns(Matrix(read_data(ticker)))
load_relative_returns_data_matrix() = get_relative_returns(Matrix(read_data()))
load_relative_centralize_returns_data_matrix() = centralize_returns(load_relative_returns_data_matrix())

# transformations
get_relative_returns(R) = diff(R,dims=1) ./ R[1:end-1,:]
get_expected_return(R)  = sum(R,dims=1)/(size(R)[1]-1)
centralize_returns(R)   = R - repeat(get_expected_return(R),size(R)[1],1)

end

## Utility functions 
# get_csv_path(str::String) = data_dir*"TIME_SERIES_INTRADAY_$(str)_60min.csv"
# get_csv_path(int::Int) = data_dir*"TIME_SERIES_INTRADAY_$(tickers[int])_60min.csv"    
# read_close_price(csv_file) =  CSV.read(csv_file, DataFrame)[:,:close] ### NOTE THAT THE TIMES MAY DIFFER


## math functions 
# """α = (α₁,α₂,...,αₙ) -> 𝔼[(r₁-μ)ᵅ¹(r₂-μ)ᵅ²⋯(rₙ-μ)ᵅⁿ] """
# function get_data_moment(R,α)
#     @assert sum(α) > 1
#     @assert sum(α) < 5
#     M = size(R)[1] 
#     α_num = findall(α .> 0)
#     return multinomial(α)*sum([prod(R[i,α_num].^ α[α_num]) for i ∈ 1:M ])/(M-1) # confusing
# end


# """α -> (∑ᵢ αᵢ )!/(α₁!⋅α₂!⋅⋅⋅αₙ!)"""
# multinomial(α) = factorial(sum(α))/prod(factorial.(α))


# """L(f(x)) = (-1)⁻¹*L(xᵀΦ⁽ᵏ⁾(x ⊗ ⋯ ⊗ x))"""
# function make_objective_function(N,k,Lx)
#     dkm  = moments.make_mon_expo(N,k,isle=false)
#     Lxk  = moments.get_Lxᵅ(Lx,dkm) 
#     R    = load_relative_returns_data_matrix([1:N...])
#     datk = map(α -> get_data_moment(R,α),dkm)
#     return (-1)^(k+1) * sum(Lxk .* datk)
# end

