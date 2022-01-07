module POP
using DataFrames

include("using_SumOfSquares.jl") ; using .using_SumOfSquares

export get_POP_estimates,
       test

function get_POP_estimates(mom_mat,heir_lvl_range)
    df = DataFrame(hier_lvl=Int[], 
               obj_val=Float64[],
               p_stat=String[],
               d_stat=String[],
               sol_time_s=Float64[])

    for s ∈ heir_lvl_range
        obj_val, p_stat, d_stat, sol_t = using_SumOfSquares.get_SOS_bound(mom_mat, s)
        push!(df,  (s, obj_val, string(p_stat), string(d_stat), sol_t))
    end
    return df
end

function test()
    error("No test coded")
end


end