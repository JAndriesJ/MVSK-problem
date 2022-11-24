module pareto_set

using DataFrames, CSV
using Dates

include("obj.jl")                   ; using .obj

export  populate_results_csv,
        read_results_csv
 
function create_empty_results_csv(save_path = pwd()*"\\"*replace(string(now()),":"=>"-")*s*".csv")
    if isfile(save_path)
        println("Using provided save path $save_path")
    else
        touch(save_path)
        println("Creating save path $save_path")
        open(save_path, "w") do io
            write(io, "λ|is_conv|w|F1|F2|F3|F4|t|stat\n")
        end;
    end
    
    return save_path
end

function populate_results_csv(R,M,V,S,K, λ_set, conv_mask; save_path="", box_size=1, sub=0)
    isfile(save_path) ? nothing : create_empty_results_csv(save_path)
    np = length(λ_set)
    F = obj.calc_F(R,M)
    for i ∈ 1:np  
        is_con = conv_mask[i]
        λ = λ_set[i]
        append_results_csv(F,M,V,S,K,λ,is_con,save_path,box_size=box_size , sub=sub)
    end
    return save_path
end
function append_results_csv(F,M,V,S,K, λ, is_con, save_path; box_size=0, sub=0)
    w, st, os, Fw = calc_results_line(F,M,V,S,K,λ, box_size=box_size, sub=sub)
    open(save_path, "a") do io
        write(io, "$(λ)|$(is_con)|$(w)|$(Fw[1])|$(Fw[2])|$(Fw[3])|$(Fw[4])|$(st)|$(os)\n")
    end
end
function calc_results_line(F,M,V,S,K,λ; box_size=0, sub=0)
    F_λ_opt_ = obj.get_F_λ_opt(M, V, S, K, λ, box_size=box_size, sub=sub)
    w  = F_λ_opt_.optimizer 
    st = F_λ_opt_.sol_time
    os = F_λ_opt_.opt_stat
    Fw =  F[1](w),F[2](w),F[3](w),F[4](w)
    return w, st, os, Fw
end

read_results_csv(csv_path) = clean_results_csv(
                                               CSV.read(csv_path, DataFrame, delim="|")
                                                )
function clean_results_csv(df)
    df.λ       = txt_to_math.(df.λ) 
    df.w       = txt_to_math.(df.w)
    return df
end
txt_to_math(s) = eval(Meta.parse(s))

end