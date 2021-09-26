module get_data
using DataFrames
using CSV 
data_dir = string(@__DIR__)*"\\assets\\"

export  download_data,
        load_data

"""Downloads time series data from "https://www.alphavantage.co/documentation/" """
function download_data()
    ts = "TIME_SERIES_INTRADAY" # TIME_SERIES_INTRADAY_EXTENDED TIME_SERIES_DAILY
    sym = "IBM" # TSCO.LON SHOP.TRT GPV.TRV DAI.DEX
    interval =  "1min"  #, 5min, 15min, 30min, 60min
    apikey = "JTJ0PUW7GZ9SAAS3"
    url = "https://www.alphavantage.co/query?function=$ts&symbol=$sym&interval=$interval&apikey=$apikey&datatype=csv"

    download(url,output=data_dir*"$(ts)_$(sym)_$(interval).csv")
end

"""Loads a predifined DataFrame"""
function load_data()
    # csv_files = readdir(data_dir,join=true)
    df = Dict()
    for sym in ["DAC" "IBM" "SAVA" "SHOP" "TSCO"]
        csv_file = "C:\\Users\\jandr\\MahCodes\\MVSK\\assets\\extended_intraday_$(sym)_15min_year1month1_adjusted.csv"
        df[sym] = CSV.read(csv_file, DataFrame)[:,[:time,:close]]
    end
    dfj = innerjoin( df["DAC"],df["IBM"],df["SAVA"],df["SHOP"],df["TSCO"], on = :time,makeunique=true)
    return rename(dfj, :close => :DAC, :close_1 =>:IBM , :close_2 => :SAVA, :close_3 => :SHOP, :close_4 => :TSCO)
end


end


