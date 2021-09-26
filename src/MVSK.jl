module MVSK

using Pkg
Pkg.status()
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Plots")
using Plots
using Statistics

using Dates

include("src\\get_data.jl")
using .get_data


df = get_data.load_data()
# Compute the empirical moments of the data (get a package that does this)
 # plot a moving window of the the moments over time
# read further on what we want to do.






using Distributions # Pkg.add("Distributions")
using StatsPlots # Pkg.add("StatsPlots")

using StatsBase ; const sb = StatsBase  # Pkg.add("StatsBase")

# All the moments are centered around their mean.
sb.moment(df.DAC, 0) 
sb.moment(df.DAC, 1) 
sb.moment(df.DAC, 2) - var(df.DAC)
sb.moment(df.DAC, 3)

end
