module stock_data
using CSV, DataFrames
data_dir = string(@__DIR__)*"\\assets\\"

export  download_data,    
        get_tickers,
        load_data,
        load_relative_returns_data_matrix,
        get_data_moment

global data_dir = pwd()*"\\assets\\"

function get_tickers()
    csv_file = data_dir*"S&P500Tickers.csv"
    return CSV.read(csv_file, DataFrame)[:,:Symbol]
end
get_tickers(range) = get_tickers()[range]
global tickers = get_tickers()



"""Loads data"""
function load_data(ticker::String)
        csv_file = get_csv_path(ticker) 
        return  read_close_price(csv_file)
end

function load_data(ticker::Int)
    csv_file = get_csv_path(tickers[ticker])
    return read_close_price(csv_file)
end

function load_data(ticker_list::Union{Vector{Int},Vector{String}})
    csv_files = [get_csv_path(tic) for tic in ticker_list]
    close_prices = [read_close_price(csv_file) for csv_file in csv_files]
    return hcat(close_prices...)
end

function load_data()
    csv_files = [get_csv_path(tic) for tic in get_tickers()]
    close_prices = [read_close_price(csv_file) for csv_file in csv_files]
    return hcat(close_prices...)
end

load_relative_returns_data_matrix(ticker) = get_relative_returns(load_data(ticker))
load_relative_returns_data_matrix() = get_relative_returns(load_data())


## Utility functions 
get_csv_path(str::String) = data_dir*"TIME_SERIES_INTRADAY_$(str)_60min.csv"
get_csv_path(int::Int) = data_dir*"TIME_SERIES_INTRADAY_$(tickers[int])_60min.csv"    
read_close_price(csv_file) =  CSV.read(csv_file, DataFrame)[:,:close] ### NOTE THAT THE TIMES MAY DIFFER

get_relative_returns(R) = diff(R,dims=1) ./ R[1:end-1,:]
get_expected_return(R) = sum(R,dims=1)/(size(R)[1]-1)
centralize_returns(R) =  R - repeat(get_expected_return(R),size(R)[1],1)

## math functions 
"""α = (α₁,α₂,...,αₙ) -> 𝔼[(r₁-μ)ᵅ¹(r₂-μ)ᵅ²⋯(rₙ-μ)ᵅⁿ] """
function get_data_moment(R,α)
    @assert sum(α) > 1
    M = size(R)[1] 
    α_num = findall(α .> 0)
    return multinomial(α)*sum([prod(R[i,α_num].^ α[α_num]) for i ∈ 1:M ])/(M-1) # confusing
end

"""α -> (∑ᵢ αᵢ )!/(α₁!⋅α₂!⋅⋅⋅αₙ!)"""
multinomial(α) = factorial(sum(α))/prod(factorial.(α))


## Getting data
"""Downloads time series data from "https://www.alphavantage.co/documentation/" 
You don't have to do this if you already have the data.
"""
function download_data()
    ticker_list = get_tickers()
    for tic in ticker_list
        ts = "TIME_SERIES_INTRADAY" # TIME_SERIES_INTRADAY_EXTENDED TIME_SERIES_DAILY
        sym = tic # TSCO.LON SHOP.TRT GPV.TRV DAI.DEX
        interval =  "60min"  #1min, 5min, 15min, 30min, 60min
        apikey = "JTJ0PUW7GZ9SAAS3"
        url = "https://www.alphavantage.co/query?function=$ts&symbol=$sym&interval=$interval&apikey=$apikey&datatype=csv"

        sleep(13) 
        download(url,data_dir*"$(ts)_$(sym)_$(interval).csv")
    end
end


end



